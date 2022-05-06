# Step 2: GSEA based on single gene survival analysis results
### Input 1: an ID conversion info file (hg19_symbol_ensembl.csv)
### Input 2: Gene set definition file in GMT format (msigdb.symbols.gmt)
### Input 3: Survival analysis result including p-values, which is output from Step 1 (Series files *_coxph_DSS_RNAseq_fpkm.csv)
### Output: GSEA analysis result for all gene sets (Series files _fgsea.csv)
rm(list = ls())

library(data.table)
library(tidyverse)
library(fgsea) ## gmtPathways, fgsea, plotEnrichment

info <- fread(
  "data/hg19_symbol_ensembl.csv", header = TRUE
) %>%
  mutate(
    ensembl = gsub("\\..*", "", ensembl)
  )

## GSEA pathway
gmt <- "data/msigdb.symbols.gmt"
pathways <- fgsea::gmtPathways(gmt)

## survival result
myfiles <- list.files(
  path = "result/",
  pattern = "*_coxph_DSS_RNAseq_fpkm.csv",
  full.names = TRUE
)

for (myfile in myfiles) {
  
  cancer <- gsub("_.*", "", basename(myfile))
  
  survres <- fread(
    myfile, data.table = FALSE
  ) %>%
    mutate(
      ensembl = gsub("\\..*", "", ensembl)
    ) %>%
    inner_join(info, by = "ensembl") %>%
    select(symbol, pvalue) %>%
    na.omit()
  
  ranks <- survres$pvalue
  names(ranks) <- survres$symbol
  
  res <- fgsea::fgsea(
    pathways = pathways,
    stats = ranks,
    nperm = 100000
  )
  
  path <- paste0("result/", cancer, "_fgsea.csv")
  fwrite(res, path)
  
  
  for (i in seq_len(nrow(res))){
    pathway <- res$pathway[[i]]
    p <- fgsea::plotEnrichment(pathways[[pathway]],ranks) + 
      labs(title=pathway)
    ggsave(paste0("result/GSEA/plotEnrichment/",pathway,"_Enrichmentplot.jpg"), p, width = 13, height = 9)
    
  }
}







