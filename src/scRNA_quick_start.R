library(DropletUtils)
library(scater)
library(tidyverse)
library(scran)
library(BiocParallel)
library(tricycle)

args = commandArgs(trailingOnly=TRUE)
hvg_num <- args[1] %>% as.integer()
pca_num <- args[2] %>% as.integer()
out_name <- args[3]

sce_full <- read10xCounts('aggr/outs/count/filtered_feature_bc_matrix.h5')
features <- read_tsv('aggr/outs/count/filtered_feature_bc_matrix/features.tsv.gz', col_names = FALSE)
mito_ens <- features %>% filter(grepl("^mt",X2)) %>% pull(X1)
is.mito <- row.names(sce_full) %in% mito_ens

qcstats <- perCellQCMetrics(sce_full, subsets=list(Mito=is.mito))

filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")


sce <- sce_full[, !filtered$discard]

# Normalization.
sce <- logNormCounts(sce)

# Feature selection.
## 10X cellranger appends a digit to the barcode to denote a different sample
colData(sce)$sample_number <- colData(sce)$Barcode  %>% str_extract('\\d+') %>% as.factor()
if (!is.na(args[4])){
	# user gives the barcode suffix(es) to filter down to
	custom_cutdown <- str_split(args[4], ',')[[1]] %>% as.integer()
	sce <- sce[,(colData(sce)$sample_number %in% custom_cutdown)]
}

## run var calc blocked on sample
dec <- modelGeneVar(sce, block = sce$sample_number, BPPARAM = MulticoreParam(4))
hvg <- getTopHVGs(dec, n = hvg_num)

# dim red
sce <- runPCA(sce, ncomponents=pca_num, subset_row=hvg)
sce <- runUMAP(sce, dimred = 'PCA')

# clustering
clusters15 <- clusterCells(sce, use.dimred="PCA", BLUSPARAM=bluster::SNNGraphParam(k=15))
clusters5 <- clusterCells(sce, use.dimred="PCA", BLUSPARAM=bluster::SNNGraphParam(k=5))
colData(sce)$cluster <- clusters15
colData(sce)$cluster5 <- clusters5

# markers
markers <- scoreMarkers(sce, groups = colData(sce)$cluster, BPPARAM = MulticoreParam(4))
markers5 <- scoreMarkers(sce, groups = colData(sce)$cluster5, BPPARAM = MulticoreParam(4))

# cell cycle
sce <- estimate_cycle_position(sce)
sce <- estimate_Schwabe_stage(sce,
                              gname.type = 'ENSEMBL',
                              species = 'mouse')
colData(sce)$tricycleStage <- colData(sce) %>%
  as_tibble() %>%
  mutate(tricycleStage = case_when(tricyclePosition < 1 ~ 'M/G1',
                                   tricyclePosition > 5 ~ 'M/G1',
                                   tricyclePosition > 2.4 & tricyclePosition < 3.2 ~ 'S',
                                   tricyclePosition > 3 & tricyclePosition < 5 ~ 'G2M',
                                   TRUE ~ 'G1/S')) %>%
  pull(tricycleStage)
# write
save(sce_full, sce, hvg, features, markers, markers5, file = out_name)
HDF5Array::saveHDF5SummarizedExperiment(x = sce, dir = gsub('.obj.Rdata', '', out_name))
