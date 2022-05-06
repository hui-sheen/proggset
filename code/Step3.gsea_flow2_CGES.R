# Step 3: Gene set survival analysis based on CGES
### Input 1: an ID conversion info file (hg19_symbol_ensembl.csv)
### Input 2: file of survival outcome (survival_clinical.csv)
### Input 3: file of gene expression (serial files *_RNASeq_fpkm.csv)
### Input 4: GSEA analysis result for all gene sets, which is output from Step 3 (Series files _fgsea.csv)
### Output:

rm(list = ls())

library(data.table)
library(tidyverse)
library(survival)
source("R/lib_proggset.R")

info <- fread(
  "data/hg19_symbol_ensembl.csv", header = TRUE
) %>%
  mutate(ensembl = gsub("\\..*", "", ensembl))

survdat <- fread(
  "data/survival_clinical.csv", data.table = FALSE
) %>%
  dplyr::select(
    patient = bcr_patient_barcode, 
    time = DSS.time,  ## OS.time, DSS.time
    status = DSS  ## OS, DSS
  ) %>%
  na.omit()

myfiles <- list.files(
  path = "data/RNA-Seq/FPKM/",
  pattern = "_RNASeq_fpkm.csv",
  full.names = TRUE
)

for (myfile in myfiles) {
  
  cancer = gsub("_.*", "", basename(myfile))
  
  fpkm <- fread(
    myfile, header = TRUE
  ) %>%
    setNames(gsub("\\..*", "", names(.)))
  
  gseafile <- paste0("result/", cancer, "_fpkm_gfsea.csv")
  if(!file.exists(gseafile)) next
  
  gseares <- fread(
    gseafile, header = TRUE
  ) %>%
    filter(
      pval < 0.05,
      NES < 0
    )
  if(nrow(gseares) == 0) next
  
  pathways <- strsplit(gseares$leadingEdge, split = "\\|")
  names(pathways) <- gseares$pathway
  
  res <- lapply(pathways, function(pathway){
    subinfo <- subset(info, symbol %in% pathway)
    
    subfpkm <- fpkm %>%
      dplyr::select(patient, any_of(subinfo$ensembl))
    
    dat <- survdat %>%
      inner_join(subfpkm, by = "patient")
    
    if(nrow(dat) == 0) return(NULL)
    
    perm_res <- calScore.perm(
      dat[,4:ncol(dat)], dat$time, dat$status, 
      quantileValues=c(0.5, 0.5), 
      nPerm = 100, seed = 100
    )
    
    r <- calScore.resub(dat[,4:ncol(dat)], dat$time, dat$status, quantileValues=c(0.5, 0.5))
   
    tryCatch(
      {
        fit <- summary(coxph(Surv(tim, cens) ~ score, data = r$res))
        
        df <- data.frame(
          cancer = cancer,
          genes = paste(subinfo$symbol[subinfo$ensembl %in% names(dat)], collapse = ","),
          beta = fit$coefficients["score", "coef"],
          HR = fit$coefficients["score", "exp(coef)"],
          lower95 = fit$conf.int["score", "lower .95"],
          upper95 = fit$conf.int["score", "upper .95"],
          pvalue = fit$coefficients["score", "Pr(>|z|)"],
          n_perm = 100,
          perm_pvalue = perm_res$perm.p.value
        )
      },
      error = function(e) return(NULL)
    )
  }) %>%
    bind_rows(.id = "pathway")
  
  path <- paste0("result/", cancer, "_GSEA_geneset_permutation_fpkm.csv")
  fwrite(res, path)
}



