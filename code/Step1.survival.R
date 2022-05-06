# Step 1: Single-gene survival analysis
### Input 1: file of survival outcome (survival_clinical.csv)
### Input 2: file of gene expression (serial files *_RNASeq_fpkm.csv)
### Output: table of Cox analysis results for all genes (serial file *_coxph_DSS_RNAseq_fpkm.csv)
rm(list = ls())

library(data.table)
library(tidyverse)
library(survival) ## Surv, coxph
library(survminer) ## ggsurvplot, surv_fit

## survival data
survdat <- fread(
  "data/survival_clinical.csv"
) %>%
  dplyr::select(
    patient = bcr_patient_barcode, 
    time = DSS.time,
    status = DSS
  ) %>%
  na.omit() %>%
  data.table::setkey("patient")

## RNAseq FPKM data
myfiles <- list.files(
  path = "data/",
  pattern = "*_RNASeq_fpkm.csv",
  full.names = TRUE
)

for (myfile in myfiles) {
  cancer <- gsub("_.*", "", basename(myfile))
  expdat <- fread(
    myfile, header = TRUE
  ) %>%
    data.table::setkey("patient")
  
  dat <- merge(survdat, expdat, by = "patient", all = FALSE)
  if(nrow(dat) == 0) next
  
  ensembls <- names(dat)[grepl("^ENSG", names(dat))]
  
  res <- lapply(ensembls, function(ensembl){
    fmla <- as.formula(paste0("Surv(time, status) ~ ", ensembl))
    
    tryCatch(
      {
        fit <- summary(survival::coxph(fmla, data = dat))
        
        ## survival plot
        dat <- dat %>%
          mutate(
            !!sym(ensembl) := ifelse(!!sym(ensembl) > median(!!sym(ensembl)), 1, 0)
          )
        fit1 <- survminer::surv_fit(fmla, data = dat)
        ggsurvplot(
          fit1, data = dat,
          title = NULL, ylab = "Disease-specific Survival Probability",
          legend.title = ensembl, legend.labs = c("low","high"),
          censor.shape = "|", censor.size = 3,
          pval = TRUE, conf.int = FALSE,
          # risk.table = FALSE, tables.height = 0.2,
          linetype = "solid",
          palette = c("#E7B800", "#2E9FDF"),
          ggtheme = theme_classic()
        )
        ggsave(filename = paste0("result/", cancer, "_", ensembl, "_DSS_survplot.jpeg"), width = 4.5, height = 3.5, dpi = 300)
        
        ## survival table
        tb <- data.frame(
          cancer = cancer,
          ensembl = ensembl,
          N = fit$n,
          events = fit$nevent,
          beta = fit$coefficients[ensembl, "coef"],
          HR = fit$coefficients[ensembl, "exp(coef)"],
          lower95 = fit$conf.int[ensembl, "lower .95"],
          upper95 = fit$conf.int[ensembl, "upper .95"],
          pvalue = fit$coefficients[ensembl, "Pr(>|z|)"]
        )
        return(tb)
      },
      error = function(e) return(NULL)
    )
  }) %>%
    bind_rows() %>%
    mutate( fdr = p.adjust(pvalue, method = "fdr") )
  
  path <- paste0("result/", cancer, "_coxph_DSS_RNAseq_fpkm.csv")
  fwrite(res, path)
}



