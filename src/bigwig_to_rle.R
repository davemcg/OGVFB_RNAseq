#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
bw_file <- args[1]
gtf_file <- args[2]
out_dir <- args[3]

library(plyranges, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(GenomicRanges, quietly = TRUE)
library(data.table, quietly = TRUE)

# 1. Setup Data & Summarize Widest Positions
# Convert to a data.frame first to safely summarize coordinates 
# and collapse duplicate gene_names before casting back to GRanges.
genes_df <- read_gff(gtf_file) %>% 
  filter(type == "gene", gene_type == 'protein_coding') %>%
  as.data.frame() %>%
  group_by(gene_name) %>%
  summarise(
    # In rare cases like PAR genes on X/Y, we take the primary chromosome listed
    seqnames = first(seqnames), 
    strand = first(strand),
    start = min(start),
    end = max(end),
    .groups = "drop"
  )

genes_gr <- makeGRangesFromDataFrame(genes_df, keep.extra.columns = TRUE)

search_ranges <- genes_gr
start(search_ranges) <- pmax(1, start(search_ranges) - 5000)
end(search_ranges) <- end(search_ranges) + 5000

sample_name <- tools::file_path_sans_ext(basename(bw_file))

# 2. Extract Data
bw_data <- read_bigwig(bw_file, overlap_ranges = search_ranges)
cov_rle_list <- coverage(bw_data, weight = "score")

# 3. Process Genes
for (i in seq_along(genes_gr)) {
  
  gene_name <- genes_gr$gene_name[i]
  chr_name <- as.character(seqnames(genes_gr[i]))
  
  extended_start <- start(search_ranges[i])
  extended_end <- end(search_ranges[i])
  
  if (chr_name %in% names(cov_rle_list)) {
    cov_rle <- cov_rle_list[[chr_name]]
    
    # Pad with zeros if chromosome falls short of extended end
    if (length(cov_rle) < extended_end) {
      cov_rle <- c(cov_rle, Rle(0, extended_end - length(cov_rle)))
    }
    
    counts <- as.numeric(cov_rle[extended_start:extended_end])
    
    # Only process and write if there is SOME signal in this range
    if (max(counts) > 0) {
      dt_long <- data.table(
        coordinate_1based = extended_start:extended_end,
        count = counts
      )
      
      # KEEP ONLY ROWS THAT ARE DIFFERENT FROM THE PREVIOUS ROW
      # fill = -1 ensures the very first row is always evaluated and kept
	  # | .I == .N ensures the LAST value is always printed
	  dt_long <- dt_long[count != shift(count, fill = -1) | .I == .N]      
      gene_dir <- file.path(out_dir, gene_name)
      
      # Safely create directory if it doesn't exist
      if (!dir.exists(gene_dir)) {
        try(dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE), silent = TRUE)
      }
      
      file_path <- file.path(gene_dir, paste0(sample_name, ".tsv.gz"))
      
      # Write sparse, compressed TSV WITHOUT headers
      fwrite(dt_long, file_path, sep = "\t", col.names = FALSE)
    }
  }
}
