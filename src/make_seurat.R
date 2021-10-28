library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(scran)
library(dyno)
library(tidyverse)
library(Seurat)
library(annotables)
library(ggridges)
library(schex)
library(clustree)
library(future)
plan(strategy = "multicore", workers = 6)
# the first term is roughly the number of MB of RAM you expect to use
# 40000 ~ 40GB
options(future.globals.maxSize = 2400000 * 1024^2)


args <- commandArgs(trailingOnly = TRUE)
rpe <- read10xCounts(args[1])

latent = args[2] %>% as.numeric()
n_neighbors = args[3] %>% as.numeric()
min_dist = args[4] %>% as.numeric()

colnames(rpe)<-colData(rpe)$Barcode # label with the barcode
row.names(rpe) <- row.names(rpe) %>% 
  enframe() %>% 
  select(ensgene = value) %>% 
  left_join(., grch38) %>%
  select(ensgene, symbol) %>% 
  unique() %>% 
  mutate(symbol = case_when(is.na(symbol) ~ ensgene, TRUE ~ symbol)) %>% 
  #mutate(out = paste(symbol, ensgene, sep='_')) %>% 
  pull(symbol) # change to gene names
row.names(rpe)[row.names(rpe) %>% duplicated()] <- paste0(row.names(rpe)[row.names(rpe) %>% duplicated()], '_1')
m <- assay(rpe,"counts") # extract counts
gg<-Matrix::rowSums(m)>0  # remove genes that are all zero
rpe<-rpe[gg,] # remove empty genes from SCE rpe 


# https://osca.bioconductor.org/workflow-integrating-datasets.html#preprocessing
rpe <- calculateQCMetrics(rpe)
# remove outliers
outlier_lib <- isOutlier(rpe$log10_total_counts, type="lower", nmad=3)
outlier_gene <- isOutlier(rpe$log10_total_features_by_counts, type="lower", nmad=3)
table(outlier_lib)
table(outlier_gene)
rpe <- rpe[,!(outlier_lib | outlier_gene)]


# normalize data
clusters = quickCluster(rpe, min.size=100)

rpe <- computeSumFactors(rpe, cluster = clusters)
rpe <- normalizeSCE(rpe)


# convert to seurat as natural log scale for counts
rpe.seurat <- as.Seurat(rpe)


# mito and type (time point) metadata
mito.genes <- grep(pattern = "^MT-", x = rownames(GetAssayData(rpe.seurat)), value = TRUE) # find mito genes
percent.mito <- Matrix::colSums(GetAssayData(rpe.seurat, slot = 'counts')[mito.genes, ])/Matrix::colSums(GetAssayData(rpe.seurat, slot = 'counts'))


rpe.seurat <- AddMetaData(object = rpe.seurat, metadata = percent.mito, col.name = "percent.mito")
type <- colnames(rpe.seurat) %>% 
                            enframe() %>% 
                            mutate(Age = case_when(grepl("1$", value) ~ 0,
                                                    grepl("2$", value) ~ 10,
                                                    grepl("3$", value) ~ 15,
                                                    grepl("4$", value) ~ 20,
                                                    grepl("5$", value) ~ 25,
                                                    grepl("6$", value) ~ 2,
                                                    grepl("7$", value) ~ 40,
                                                    grepl("8$", value) ~ 60,
                                                    TRUE ~ 1000)) %>% 
                            pull(Age)
names(type) <-colnames(rpe.seurat)
rpe.seurat <- AddMetaData(object = rpe.seurat, 
                          metadata = type, 
                          col.name = 'Age')


rpe.seurat.preFilter <- rpe.seurat
# remove cells that are > 5% mito
rpe.seurat <- subset(rpe.seurat, subset = percent.mito < 0.05)
# overall remove all cells < 500 nCounts
rpe.seurat <- subset(rpe.seurat, subset = nFeature_RNA > 750 )

# dim reduction
rpe.seurat <- FindVariableFeatures(rpe.seurat)
rpe.seurat <- ScaleData(object = rpe.seurat, vars.to.regress = c("nCount_RNA", "percent.mito"))
rpe.seurat <- RunPCA(object = rpe.seurat, pc.genes = rpe.seurat@var.genes, do.print = TRUE, pcs.print = 1:5, 
                     genes.print = 5, dims = 100)

#rpe.seurat <- RunTSNE(rpe.seurat, reduction.use = "pca", dims.use = 1:20,  do.fast = TRUE, do.label =T)

rpe.seurat <- RunUMAP(rpe.seurat, reduction = "pca", dims = 1:latent, min.dist = min_dist, n.neighbors = n_neighbors)
rpe.seurat <- RunUMAP(rpe.seurat, reduction = "pca", dims = 1:latent, min.dist = min_dist, n.neighbors = n_neighbors, n.components = 3, reduction.name = 'UMAP3D', reduction.key = 'UMAP3D_')

# clustering 
rpe.seurat <- FindNeighbors(rpe.seurat, reduction = 'pca', dims = 1:latent)
rpe.seurat <- FindClusters(rpe.seurat)

save(rpe.seurat, rpe.seurat.preFilter, file = args[5], compress = FALSE)
