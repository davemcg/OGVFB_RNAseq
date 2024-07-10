library(DropletUtils)
library(scater)
library(tidyverse)
library(scran)
library(BiocParallel)
library(tricycle)

args = commandArgs(trailingOnly=TRUE)
ssheet <- read_csv(args[1])
hvg_num <- args[2] %>% as.integer()
pca_num <- args[3] %>% as.integer()
out_name <- args[4]

sce_full <- read10xCounts('aggr/outs/count/filtered_feature_bc_matrix.h5')
features <- read_tsv('aggr/outs/count/filtered_feature_bc_matrix/features.tsv.gz', col_names = FALSE)
features$X2 <- toupper(features$X2)
mito_ens <- features %>% filter(grepl("^MT",X2)) %>% pull(X1)
is.mito <- row.names(sce_full) %in% mito_ens

solo_scores <- read_csv('aggr/solo.csv.gz')

colData(sce_full)$sample_number <- colData(sce_full)$Barcode  %>% str_extract('\\d+') %>% as.factor()


# loop through each individual sample to identify outliers with the scuttle toolset
filter_list <- list()
for (i in ssheet$sample_id){
	sce_sub <-  read10xCounts(paste0(i,'/outs/filtered_feature_bc_matrix.h5'))
	qcstats_sub <- perCellQCMetrics(sce_sub, subsets=list(Mito=is.mito))
	filter_list[[i]] <- quickPerCellQC(qcstats_sub, percent_subsets="subsets_Mito_percent") %>% data.frame()
}
filtered <- data.table::rbindlist(filter_list)

qcstats <- perCellQCMetrics(sce_full, subsets=list(Mito=is.mito))
sce_full <- addPerCellQCMetrics(sce_full, subsets=list(Mito=is.mito))
sce <- sce_full[, !filtered$discard]

# Normalization.
if (!is.na(args[5])){
	if (grepl("downsample", args[5])){
		sce <- logNormCounts(sce, downsample = TRUE)
	} else {
		sce <- logNormCounts(sce)
	}
} else {
	sce <- logNormCounts(sce)
}

# Feature selection.
## 10X cellranger appends a digit to the barcode to denote a different sample
colData(sce)$sample_number <- colData(sce)$Barcode  %>% str_extract('\\d+') %>% as.factor()
if (!is.na(args[4])){
	# user gives the barcode suffix(es) to filter down to
	custom_cutdown <- str_split(args[5], ',')[[1]] %>% as.integer()
	sce <- sce[,(colData(sce)$sample_number %in% custom_cutdown)]
}

## run var calc blocked on sample
dec <- modelGeneVar(sce, block = sce$sample_number, BPPARAM = MulticoreParam(4))
hvg <- getTopHVGs(dec, n = hvg_num)

# remove mito or ribo genes from hvg
hvg_sub <- features %>% filter(X1 %in% hvg) %>% filter(!grepl("^MT|^RPS|^RPL", X2)) %>% pull(X1)
# dim red
sce <- runPCA(sce, ncomponents=pca_num, subset_row=hvg_sub, exprs_values = "logcounts", scale = TRUE)
sce <- runUMAP(sce, dimred = 'PCA')

# clustering
print("clustering")
clusters15 <- clusterCells(sce, use.dimred="PCA", BLUSPARAM=bluster::SNNGraphParam(k=15))
clusters5 <- clusterCells(sce, use.dimred="PCA", BLUSPARAM=bluster::SNNGraphParam(k=5))
colData(sce)$cluster <- clusters15
colData(sce)$cluster5 <- clusters5

# markers
print("markers")
markers <- scoreMarkers(sce, groups = colData(sce)$cluster, BPPARAM = MulticoreParam(4))
markers5 <- scoreMarkers(sce, groups = colData(sce)$cluster5, BPPARAM = MulticoreParam(4))

# cell cycle
print("cell cycle")
guess_species <- ifelse(grepl("ENSG", features$X1[1]), "human", "mouse")
print("estimate_cycle_position")
sce <- try({estimate_cycle_position(sce,
                              gname.type = 'ENSEMBL',
                              species = guess_species) })
print("get stage")
colData(sce)$tricycleStage <- try({ colData(sce) %>%
  as_tibble() %>%
  mutate(tricycleStage = case_when(tricyclePosition < 1 ~ 'M/G1',
                                   tricyclePosition > 5 ~ 'M/G1',
                                   tricyclePosition > 2.4 & tricyclePosition < 3.2 ~ 'S',
                                   tricyclePosition > 3 & tricyclePosition < 5 ~ 'G2M',
                                   TRUE ~ 'G1/S')) %>%
  pull(tricycleStage) })
# write
save(sce_full, qcstats, solo_scores, filtered, sce, hvg, hvg_sub, features, markers, markers5, file = out_name)
HDF5Array::saveHDF5SummarizedExperiment(x = sce, dir = gsub('.obj.Rdata', '', out_name))
