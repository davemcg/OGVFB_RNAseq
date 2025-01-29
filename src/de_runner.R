de_runner <- function(txi,
                      metadata,
                      design = "~batch + treatment",
                      reduced = "~batch",
                      ntop = 1000,
                      lrt = TRUE,
                      ensembl = TRUE,
                      deseq2_input = TRUE){
  if (!deseq2_input){
    dds_temp <- DESeqDataSetFromTximport(txi,
                                         colData = metadata,
                                         design = as.formula(design))}
  else {dds_temp <- txi}

  if (lrt){
    dds_temp <- DESeq(dds_temp, test = 'LRT',
                      reduced = as.formula(reduced),
                      parallel = TRUE)
  } else {
    dds_temp <- DESeq(dds_temp,
                      parallel = TRUE)
  }

  if (ensembl){
    res_temp <- results(dds_temp) %>%
      as_tibble(rownames = 'ENSEMBL') %>%
      arrange(padj) %>%
      left_join(conv_table) %>%
      relocate(SYMBOL)
  } else{
    res_temp <- results(dds_temp) %>%
      as_tibble(rownames = 'SYMBOL') %>%
      arrange(padj)
  }
  out <- list()
  out$dds <- dds_temp
  out$res <- res_temp
  # corrected counts
  vst<-vst(dds_temp,blind = FALSE)
  #assay(vst) <- limma::removeBatchEffect(assay(vst), batch = vst$batch, batch2 = vst$leiden2)
  out$vst <- assay(vst)
  out$counts <- assay(dds_temp, 'counts')
  out$cpm <- metamoRph::normalize_data(assay(dds_temp, 'counts'))
  ## pca ##
  ntop = ntop
  Pvars <- rowVars(out$vst)
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,
                                                        length(Pvars)))]
  PCA <- prcomp(t(out$vst[select, ]), scale = F)
  PCA$percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  out$PCA <- PCA

  out
}

go_runner <- function(res,
                      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                      orgdb_chr = "org.Hs.eg.db",
                      fromType="SYMBOL",
                      toType="ENTREZID",
                      res_padj_threshold = 0.01,
                      logFC_threshold = 1,
                      mode = 'GO',
                      ...){

  diff_genes <- res %>%
    filter(padj < res_padj_threshold,
           !grepl("RPS|RPL",Gene),
           abs(log2FoldChange) > logFC_threshold)

  eg_diff_genes <- bitr(diff_genes$Gene, fromType=fromType, toType=toType, OrgDb="org.Hs.eg.db")
  eg_diff_genes <- diff_genes %>% left_join(., eg_diff_genes, by = c('Gene'=fromType))

  eg_universe = bitr(res %>% filter(!grepl("RPS|RPL",Gene)) %>%  pull(Gene), fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb_chr)
  if (mode == 'GO'){
    enrichGO(gene          = eg_diff_genes$ENTREZID,
             universe      = eg_universe$ENTREZID,
             OrgDb         = OrgDb,
             ont           = "all",
             pvalueCutoff = Inf,
             qvalueCutoff = Inf,
             pAdjustMethod = "BH",
             readable      = TRUE)
  } else if (mode == 'pathway'){
    enrichPathway(gene          = eg_diff_genes$ENTREZID,
                  universe      = eg_universe$ENTREZID,
                  pvalueCutoff = Inf,
                  qvalueCutoff = Inf,
                  pAdjustMethod = "BH",
                  readable      = TRUE)
  }
}

gsea_runner <- function(res,
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                        mode = 'GSEA',
                        ...){
  if (mode == 'GSEA'){
    gse <- gseGO(geneList=logFC_maker(res),
                 ont ="ALL",
                 keyType = "ENTREZID", pvalueCutoff = Inf,
                 OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                 pAdjustMethod = "BH",
                 eps = 0)
  } else if (mode == 'pathway'){
    gse <- gsePathway(geneList=logFC_maker(res),
                      pvalueCutoff = Inf,
                      organism ='human',
                      pAdjustMethod = "BH",
                      eps = 0)
  }
  gse <- setReadable(gse, OrgDb = OrgDb)
  gse
}

logFC_maker <- function(res,
                        fromType="SYMBOL",
                        toType="ENTREZID",
                        orgdb_chr = "org.Hs.eg.db"){
  all_genes <- bitr(res %>% filter(!grepl("RPS|RPL",Gene)) %>% pull(Gene), fromType=fromType, toType=toType, OrgDb=orgdb_chr)
  all_genes <- all_genes %>% left_join(res, by = c('SYMBOL' = 'Gene'))
  logFC <- all_genes$log2FoldChange
  names(logFC) <- all_genes$ENTREZID
  logFC <- na.omit(logFC)
  logFC = sort(logFC, decreasing = TRUE)
  logFC
}
