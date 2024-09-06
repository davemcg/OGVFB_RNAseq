# lightly adapted from ~/git/EiaD_build/src/make_counts.R
print("go")
library(tidyverse)
library(tximport)
library(rjson)
library(data.table)
print("library load")
args <- commandArgs(trailingOnly=TRUE)

files <- args
print("load files")
# process metadata
log_files <- gsub('quant.sf', 'aux_info/meta_info.json', files)

log_list <- map(log_files, function(x) fromJSON(file = x))
pos_processed <- grep('num_processed', names(log_list[[1]]) )
pos_mapped <- grep('num_mapped', names(log_list[[1]]) )
pos_date <- grep('end_time', names(log_list[[1]]) )
pos_version <- grep('salmon_version', names(log_list[[1]]) )
run_meta <- cbind(log_files,
				  map(log_list, pos_processed) %>% unlist(),
				  map(log_list, pos_mapped) %>% unlist(),
				  map(log_list, pos_date) %>% unlist(),
                  map(log_list, pos_version) %>% unlist()) %>% as_tibble() 
colnames(run_meta) <- c("log", "processed", "mapped", "date", "version") 
run_meta <- run_meta %>%
			mutate(processed = as.integer(processed),
					mapped = as.integer(mapped),
					align_perc = mapped/processed) %>%
			relocate(log, align_perc)
print('two')
# pull anno
salmon_anno  <- fread(files[1])
salmon_anno2 <- salmon_anno %>% 
					separate(Name,into = c('transcript_id','gene_id','havana_gene',
											'havana_transcript','transcript_name',
											'gene_name','length','type','empty'), 
								sep = '\\|', remove = FALSE) 

salmon_tx2gene <- salmon_anno2 %>% dplyr::select(Name, gene_id)

# import
txi <- tximport(files, type = "salmon", tx2gene = salmon_tx2gene)
txi_tx <- tximport(files, type = "salmon", tx2gene = salmon_tx2gene, txOut = TRUE)
# extract counts
puller <- function(txi_object, slot = 'counts', ncolnames = gsub('salmon_quant|\\/|quant\\.sf', '', files)){
	out <- txi_object[[slot]]
	colnames(out) <- ncolnames
	out
}

gene_counts <- puller(txi, 'counts')
gene_TPM <- puller(txi, 'abundance')

tx_counts <- puller(txi_tx, 'counts')
tx_TPM <- puller(txi_tx, 'abundance')

system('mkdir -p counts')
write_csv(salmon_anno2 %>% as_tibble(), file = 'counts/gene_anno.csv.gz')
write_csv(gene_counts %>% as_tibble(rownames = 'Gene'), file = 'counts/gene_counts.csv.gz')
write_csv(gene_TPM %>% as_tibble(rownames = 'Gene'), file = 'counts/gene_tpm.csv.gz')
write_csv(tx_counts %>% as_tibble(rownames = 'Transcript'), file = 'counts/tx_counts.csv.gz')
write_csv(tx_TPM %>% as_tibble(rownames = 'Transcript'), file = 'counts/tx_tpm.csv.gz')

write_csv(run_meta, file = 'run_meta.csv.gz')

