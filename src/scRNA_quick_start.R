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
cellbender_input <- args[5]
subsets <- args[6]


if (toupper(cellbender_input) == 'TRUE'){
  ssheet$molecule_h5 <- sub("molecule_info.h5", "cellbender_filtered.h5", ssheet$molecule_h5)
  ssheet$sample_number <- seq(1,nrow(ssheet))
  cellbender_files <- gsub("molecule_info.h5", "cellbender_filtered.h5", ssheet$molecule_h5)
  sce_list <- list()
  counter = 1
  for (i in cellbender_files){
    sce_list[[i]] <- read10xCounts(i)
    colData(sce_list[[i]])$sample_number <- counter
    counter <- counter + 1
  }
  sce_full <- Reduce(cbind, sce_list)
} else {
  sce_full <- read10xCounts('aggr/outs/count/filtered_feature_bc_matrix.h5')
  ssheet$molecule_h5 <- sub("molecule_info.h5", "filtered_feature_bc_matrix.h5", ssheet$molecule_h5)
  # Feature selection.
  ## 10X cellranger appends a digit to the barcode to denote a different sample
  colData(sce_full)$sample_number <- colData(sce_full)$Barcode  %>% str_extract('\\d+') %>% as.factor()
}

features <- read_tsv('aggr/outs/count/filtered_feature_bc_matrix/features.tsv.gz', col_names = FALSE)
#features$X2 <- toupper(features$X2)

mito_ens <- features %>% mutate(X2 = toupper(X2)) %>% filter(grepl("^MT",X2)) %>% pull(X1)
is.mito <- row.names(sce_full) %in% mito_ens

solo_scores <- read_csv('aggr/solo.csv.gz')

colData(sce_full)$sample_number <- colData(sce_full)$Barcode  %>% str_extract('\\d+') %>% as.factor()


# loop through each individual sample to identify outliers with the scuttle toolset
filter_list <- list()
for (i in seq(1,nrow(ssheet))){
  id = ssheet[i,] %>% pull(1)
  h5 <- ssheet[i,] %>% pull(2)
  sce_sub <-  read10xCounts(h5)
  qcstats_sub <- perCellQCMetrics(sce_sub, subsets=list(Mito=is.mito))
  filter_list[[id]] <- quickPerCellQC(qcstats_sub, percent_subsets="subsets_Mito_percent") %>% data.frame()
}
filtered <- data.table::rbindlist(filter_list)

qcstats <- perCellQCMetrics(sce_full, subsets=list(Mito=is.mito))
sce_full <- addPerCellQCMetrics(sce_full, subsets=list(Mito=is.mito))
sce <- sce_full[, !filtered$discard]




run_workflow <- function(sce){
  # Normalization.
  if (!is.na(subsets)){
    if (grepl("downsample",subsets)){
      sce <- logNormCounts(sce, downsample = TRUE)
    } else {
      sce <- logNormCounts(sce)
    }
  } else {
    sce <- logNormCounts(sce)
  }


  if (!is.na(out_name)){
    # user gives the barcode suffix(es) to filter down to
    custom_cutdown <- str_split(subsets, '-')[[1]] %>% as.integer()
    sce <- sce[,(colData(sce)$sample_number %in% custom_cutdown)]
  }

  # Feature selection.
  ## 10X cellranger appends a digit to the barcode to denote a different sample
  colData(sce)$sample_number <- colData(sce)$Barcode  %>% str_extract('\\d+') %>% as.factor()
  if (!is.na(out_name)){
    # user gives the barcode suffix(es) to filter down to
    custom_cutdown <- str_split(subsets, ',')[[1]] %>% as.integer()
    sce <- sce[,(colData(sce)$sample_number %in% custom_cutdown)]
  }

  ## run var calc blocked on sample
  dec <- modelGeneVar(sce, block = sce$sample_number, BPPARAM = MulticoreParam(4))
  hvg <- getTopHVGs(dec, n = hvg_num)

  # remove mito or ribo genes from hvg
  hvg_sub <- features %>% filter(X1 %in% hvg) %>% mutate(X2 = toupper(X2)) %>% filter(!grepl("^MT|^RPS|^RPL", X2)) %>% pull(X1)
  # dim red
  sce <- runPCA(sce, ncomponents=pca_num, subset_row=hvg_sub, exprs_values = "logcounts", scale = TRUE)
  sce <- runUMAP(sce, dimred = 'PCA', min_dist = 0.3, metric = 'cosine')

  # clustering
  print("clustering")
  clusters15 <- clusterCells(sce, use.dimred="PCA", BLUSPARAM=bluster::SNNGraphParam(k=15, cluster.fun = 'leiden'))
  clusters5 <- clusterCells(sce, use.dimred="PCA", BLUSPARAM=bluster::SNNGraphParam(k=5, cluster.fun = 'leiden'))
  colData(sce)$cluster <- clusters15
  colData(sce)$cluster5 <- clusters5

  # markers
  print("markers")
  markers <- scoreMarkers(sce, groups = colData(sce)$cluster, BPPARAM = MulticoreParam(4))
  #markers5 <- scoreMarkers(sce, groups = colData(sce)$cluster5, BPPARAM = MulticoreParam(4))

  fmarkers <- scran::findMarkers(sce, groups = colData(sce)$cluster, BPPARAM = MulticoreParam(4))

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

  out <- list()
  out$sce <- sce
  out$markers <- markers
  out$fmarkers <- fmarkers
  out$hvg_sub <- hvg_sub
  out$hvg <- hvg
  out
}

# run
sce <- run_workflow(sce)
sce_full <- run_workflow(sce_full)

# write
save(qcstats, solo_scores, filtered, sce, sce_full, file = out_name)
HDF5Array::saveHDF5SummarizedExperiment(x = sce$sce, dir = gsub('.obj.Rdata', '', out_name))
HDF5Array::saveHDF5SummarizedExperiment(x = sce_full$sce, dir = paste0(gsub('.obj.Rdata', '', out_name), "_all"))
