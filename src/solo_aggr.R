library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

ssheet <- read_csv(args[1]) %>% mutate(num = row_number())

solo_list <- list()

for (index in ssheet$num){
	index = as.character(index)
	sid = ssheet %>% filter(num == index) %>% pull(sample_id)
	snum = ssheet %>% filter(num == index) %>% pull(num)
	solo_list[[sid]] <- read_csv(paste0(sid,'/outs/solo.csv.gz')) %>% 
					mutate(bc = str_extract(`...1`, '[ACGT]+'), 
						num = snum, 
						bc2 = paste0(bc,'-',num))
}

solo_list %>% bind_rows() %>% dplyr::select(barcode = bc2, solo_doublet, solo_score) %>%
	write_csv(file = args[2])


