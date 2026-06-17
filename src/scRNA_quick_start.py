#!/usr/bin/env python
"""
GPU single-cell processing pipeline (rapids-singlecell port of scRNA_quick_start.R).

Usage:
    python scRNA_quick_start_gpu.py <samplesheet.csv> <hvg_num> <pca_num> \
        <out_name.h5ad> <cellbender_input(TRUE/FALSE)> <subsets> [n_gpu_streams]

Notes:
  - Targets rapids-singlecell (rsc) with a scanpy CPU fallback for steps rsc lacks.
  - Input branch implemented: CellBender per-sample .h5 (cellbender_input=TRUE).
  - Clustering: Leiden (matches R original; SNN graph in R, kNN graph here).
  - QC: per-sample adaptive (MAD) reimplementing scuttle::quickPerCellQC, OR'd with hard cutoffs.
"""

import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import rapids_singlecell as rsc
import cupy as cp
from scipy import sparse

# ----------------------------------------------------------------------
# Args
# ----------------------------------------------------------------------
ssheet_path = sys.argv[1]
hvg_num     = int(sys.argv[2])
pca_num     = int(sys.argv[3])
out_name    = sys.argv[4]
cellbender_input = sys.argv[5].upper() == "TRUE"
subsets     = sys.argv[6] if len(sys.argv) > 6 and sys.argv[6] not in ("NA", "None", "") else None

assert hvg_num > 0 and pca_num > 0, "hvg_num and pca_num must be positive integers"

ssheet = pd.read_csv(ssheet_path)
RUN_FULL = False   # mirror the optional pre-QC full run from the R script

# ----------------------------------------------------------------------
# Read CellBender per-sample .h5 and concatenate
# ----------------------------------------------------------------------
def read_cellbender_samples(ssheet):
    # column 1 = sample id, column 2 = path to molecule_info.h5 (rewritten to cellbender output)
    ssheet = ssheet.copy()
    mol_col = ssheet.columns[1]
    id_col  = ssheet.columns[0]
    ssheet[mol_col] = ssheet[mol_col].str.replace(
        "molecule_info.h5", "cellbender_filtered.h5", regex=False
    )
    adatas = []
    for i, row in ssheet.reset_index(drop=True).iterrows():
        sample_number = i + 1                      # 1-based, matches R loop index
        a = sc.read_10x_h5(row[mol_col])
        a.var_names_make_unique()
        a.obs["sample_number"] = sample_number
        a.obs["sample_id"] = str(row[id_col])
        # unique barcodes across samples: barcode + sample suffix
        a.obs_names = [f"{bc}-{sample_number}" for bc in a.obs_names]
        adatas.append(a)
    full = ad.concat(adatas, join="outer", label=None, index_unique=None)
    full.obs_names_make_unique()
    return full

if cellbender_input:
    adata_full = read_cellbender_samples(ssheet)
else:
    raise NotImplementedError("Only the CellBender branch was requested/implemented.")

# Ensure raw integer counts are preserved before any normalization
adata_full.layers["counts"] = adata_full.X.copy()

# CuPy sparse matrices reject integer dtypes (only bool/float32/float64/complex).
# 10x/CellBender .h5 counts load as int -> cast X and the counts layer to float32 before GPU transfer.
adata_full.X = adata_full.X.astype(np.float32)
adata_full.layers["counts"] = adata_full.layers["counts"].astype(np.float32)

# ----------------------------------------------------------------------
# Mito flag (gene symbols start with MT, case-insensitive)
# ----------------------------------------------------------------------
# 10x var typically has 'gene_ids' (ENSEMBL) + var_names (symbols). Guard both.
symbols = adata_full.var_names.astype(str)
adata_full.var["mt"] = symbols.str.upper().str.startswith("MT")
adata_full.var["ribo"] = symbols.str.upper().str.match(r"^(RPS|RPL)")

# ----------------------------------------------------------------------
# Per-cell QC metrics (computed on GPU)
# ----------------------------------------------------------------------
rsc.get.anndata_to_GPU(adata_full)
rsc.pp.calculate_qc_metrics(adata_full, qc_vars=["mt"])
# rsc stores: n_genes_by_counts, total_counts, pct_counts_mt
rsc.get.anndata_to_CPU(adata_full)

# ----------------------------------------------------------------------
# Per-sample adaptive QC  (reimplements scuttle::quickPerCellQC)
#   - low total_counts:      median - 3*MAD on log scale
#   - low n_genes_by_counts: median - 3*MAD on log scale
#   - high pct_counts_mt:    median + 3*MAD on natural scale
# ----------------------------------------------------------------------
def is_outlier_mad(values, nmads=3.0, log=False, lower=True, upper=True):
    v = np.asarray(values, dtype=float)
    x = np.log1p(v) if log else v
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    mad = mad if mad > 0 else 1e-9
    lo = med - nmads * mad
    hi = med + nmads * mad
    out = np.zeros(len(x), dtype=bool)
    if lower:
        out |= x < lo
    if upper:
        out |= x > hi
    return out

obs = adata_full.obs
discard = np.zeros(adata_full.n_obs, dtype=bool)
for s in obs["sample_number"].unique():
    m = (obs["sample_number"] == s).values
    sub = obs.loc[m]
    low_counts = is_outlier_mad(sub["total_counts"],          log=True,  lower=True,  upper=False)
    low_genes  = is_outlier_mad(sub["n_genes_by_counts"],     log=True,  lower=True,  upper=False)
    high_mito  = is_outlier_mad(sub["pct_counts_mt"],         log=False, lower=False, upper=True)
    discard[m] = low_counts | low_genes | high_mito

# Hard cutoffs (matches R: sum < 500 OR detected < 200)
discard |= (obs["total_counts"].values < 500) | (obs["n_genes_by_counts"].values < 200)

adata_full.obs["discard"] = discard
print(f"Discarding {discard.sum()} / {adata_full.n_obs} cells "
      f"({100*discard.mean():.1f}%) across {obs['sample_number'].nunique()} samples")

adata = adata_full[~discard].copy()

# ----------------------------------------------------------------------
# Workflow
# ----------------------------------------------------------------------
def run_workflow(a, hvg_num, pca_num, subsets):
    a = a.copy()

    # optional subset to given sample numbers e.g. "1-2-3"
    if subsets is not None and "downsample" not in subsets:
        keep = [int(x) for x in subsets.split("-")]
        a = a[a.obs["sample_number"].isin(keep)].copy()

    rsc.get.anndata_to_GPU(a, convert_all=True)   # convert_all moves .X AND .layers (incl. 'counts') to GPU; seurat_v3 HVG reads the counts layer

    # Normalization (log1p of CPM-like). R used logNormCounts (size-factor); this is the scanpy analogue.
    rsc.pp.normalize_total(a, target_sum=1e4)
    rsc.pp.log1p(a)

    # HVGs blocked by sample. seurat_v3 expects raw counts, so point it at the counts layer.
    rsc.pp.highly_variable_genes(
        a, n_top_genes=hvg_num, flavor="seurat_v3",
        batch_key="sample_number", layer="counts",
    )
    # Exclude mito/ribo from the HVG set fed to PCA (matches R behavior, done pre-PCA)
    a.var.loc[a.var["mt"] | a.var["ribo"], "highly_variable"] = False

    # PCA on HVGs, scaled (R used scale=TRUE on logcounts subset to HVGs)
    rsc.pp.scale(a, max_value=10, mask_obs=None)
    rsc.pp.pca(a, n_comps=pca_num, mask_var="highly_variable")

    # Neighbors + UMAP (R: cosine metric, min_dist 0.3)
    rsc.pp.neighbors(a, n_neighbors=15, use_rep="X_pca", metric="cosine")
    rsc.tl.umap(a, min_dist=0.3)

    # Clustering: Leiden (matches R original, which used Leiden on an SNN graph)
    print("clustering")
    rsc.tl.leiden(a)
    a.obs["cluster"] = a.obs["leiden"]

    # report cluster count and cells per cluster
    counts = a.obs["cluster"].value_counts().sort_index()
    print(f"{counts.shape[0]} clusters; cells per cluster:")
    print(counts.to_string())

    rsc.get.anndata_to_CPU(a, convert_all=True)   # bring .X and all layers back to CPU; rank_genes_groups stores results in .uns either way

    # Markers (one-vs-rest): standard rank_genes_groups, GPU Wilcoxon when available.
    # Stored under a dedicated .uns key so the later pairwise loop doesn't overwrite it.
    print("markers (one-vs-rest)")
    _rank_genes_groups(a, groupby="cluster", method="wilcoxon", tie_correct=True,
                       key_added="rank_genes_groups")
    diff = sc.get.rank_genes_groups_df(a, group=None, key="rank_genes_groups")
    if "gene_ids" in a.var.columns:
        ens_map = dict(zip(a.var_names, a.var["gene_ids"]))
        diff["ENSEMBL"] = diff["names"].map(ens_map)
    diff = diff.rename(columns={"names": "SYMBOL"})

    # Markers (all unordered pairs): cluster-vs-cluster contrasts.
    print("markers (all pairwise)")
    diff_pairwise = pairwise_rank_genes(a, groupby="cluster", method="wilcoxon")
    if "gene_ids" in a.var.columns:
        diff_pairwise["ENSEMBL"] = diff_pairwise["names"].map(ens_map)
    diff_pairwise = diff_pairwise.rename(columns={"names": "SYMBOL"})

    # Cell cycle:
    #   NOTE: this is NOT tricycle. rapids/scanpy has no tricycle equivalent.
    #   Using Scanpy's score_genes_cell_cycle (Regev S/G2M lists) -> discrete phase.
    print("cell cycle")
    s_genes, g2m_genes = _cc_gene_lists(a)
    if s_genes and g2m_genes:
        try:
            sc.tl.score_genes_cell_cycle(a, s_genes=s_genes, g2m_genes=g2m_genes)
            # a.obs['phase'] in {'G1','S','G2M'}
        except Exception as e:
            print(f"cell cycle scoring failed, skipping: {e}")

    return {"adata": a, "diff": diff, "diff_pairwise": diff_pairwise}

def _rank_genes_groups(a, **kwargs):
    """Call rsc.tl.rank_genes_groups (GPU) if available in this build, else scanpy CPU."""
    fn = getattr(getattr(rsc, "tl", None), "rank_genes_groups", None)
    if fn is not None:
        # rsc requires data on GPU; transfer .X only for the test
        rsc.get.anndata_to_GPU(a)
        fn(a, **kwargs)
        rsc.get.anndata_to_CPU(a)
    else:
        sc.tl.rank_genes_groups(a, **kwargs)


def pairwise_rank_genes(a, groupby="cluster", method="wilcoxon", tie_correct=True):
    """All unordered-pairs pairwise differential expression.

    Loops over each unordered cluster pair {A, B} once (A before B in category order)
    and runs rank_genes_groups with groups=[A], reference=B, i.e. 'cluster A vs cluster B'.
    Results are stacked long-form with explicit 'cluster' (=A, the tested group) and
    'reference' (=B) columns. logfoldchanges are A relative to B; for the B-vs-A view,
    negate logfoldchanges and swap the labels.

    For k clusters this is k*(k-1)/2 test calls. To avoid paying a GPU<->CPU transfer
    per call, the object is moved to GPU ONCE before the loop and back ONCE after when
    the GPU rank_genes_groups is available; otherwise it stays on CPU (scanpy).
    """
    cats = list(a.obs[groupby].astype("category").cat.categories)

    gpu_fn = getattr(getattr(rsc, "tl", None), "rank_genes_groups", None)
    use_gpu = gpu_fn is not None
    if use_gpu:
        rsc.get.anndata_to_GPU(a)          # one transfer for the whole loop

    frames = []
    try:
        for i in range(len(cats)):
            for j in range(i + 1, len(cats)):     # unordered: only i < j
                A, B = cats[i], cats[j]
                kw = dict(groupby=groupby, groups=[A], reference=B, method=method,
                          key_added="rank_genes_groups_pairwise")
                if method == "wilcoxon":
                    kw["tie_correct"] = tie_correct
                if use_gpu:
                    gpu_fn(a, **kw)            # already on GPU, no per-call transfer
                else:
                    sc.tl.rank_genes_groups(a, **kw)
                df = sc.get.rank_genes_groups_df(a, group=A, key="rank_genes_groups_pairwise")
                df = df.head(100)          # top 100 genes per pair (already ranked by the test score)
                df["cluster"] = A
                df["reference"] = B
                frames.append(df)
    finally:
        if use_gpu:
            rsc.get.anndata_to_CPU(a)         # single transfer back, even if the loop errors

    out = pd.concat(frames, ignore_index=True)
    return out


def _cc_gene_lists(a):
    """Regev lab S and G2M gene symbols, intersected with this object's var_names."""
    s_genes = ["MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
               "CDCA7","DTL","PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP","RAD51AP1",
               "GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51",
               "RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1",
               "CLSPN","POLA1","CHAF1B","BRIP1","E2F8"]
    g2m_genes = ["HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
                 "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2",
                 "CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1",
                 "KIF20B","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1",
                 "NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1",
                 "ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"]
    present = set(a.var_names.str.upper())
    sym_upper = {v.upper(): v for v in a.var_names}
    s = [sym_upper[g] for g in s_genes if g in present]
    g = [sym_upper[g] for g in g2m_genes if g in present]
    return s, g

# ----------------------------------------------------------------------
# Run + write
# ----------------------------------------------------------------------
res = run_workflow(adata, hvg_num, pca_num, subsets)

res["adata"].write_h5ad(out_name, compression="gzip")
res["diff"].to_csv(out_name.replace(".h5ad", ".markers.csv"), index=False)
res["diff_pairwise"].to_csv(out_name.replace(".h5ad", ".markers_pairwise.csv"), index=False)

if RUN_FULL:
    res_full = run_workflow(adata_full, hvg_num, pca_num, subsets)
    res_full["adata"].write_h5ad(out_name.replace(".h5ad", "_all.h5ad"), compression="gzip")
