library(Seurat)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

load(args[1])
dims = args[2] %>% as.numeric()
nneighbors = args[3] %>% as.numeric()
dist = args[4] %>% as.numeric()

orig_meta <- rpe.seurat@meta.data 
umap <- Embeddings(rpe.seurat[['umap']]) %>% as_tibble(rownames = 'Barcode') %>% 
  left_join(., orig_meta) %>% 
  mutate(Dims = dims,
         Neighbors = nneighbors,
         Dist = dist)
save(umap, file = args[5])
