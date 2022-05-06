#Step4
## A - DCGL + CSPN
### INPUT 1: FPKM for multiple cancer types
### INPUT 2: Sample groupping of multiple cancer types
### INPUT 3: Pathway definition identified with symbols
### OUTPUT 1: A background gene-gene network for CSPN use (DCLi.all.txt & DCLd.all.txt)
################ This file goes to cspn/reference/PPI_network
### OUTPUT 2: Gene set definition identified with EntrezGeneIDs for CSPN use (genes_hsa####.txt files)
################ These files go to cspn/reference/original_gene_sets/
### OUTPUT 3: List of concerned pathways (ccPWs.txt)
################ This file goes to cspn/reference/concerned_pathway_list/ 
#setwd('D:/WORK/crosstalk/memberTalk')
library(DCGL)
library(org.Hs.eg.db)
source("R/lib_proggset.R")
load('DCEAinit.RData') # Needs to reduce to minimal dat {fpkm}
load('syms_in_pw.RData') # list of involved genes; Needs to reduce to entailed pathways
kept_sym <- intersect(rownames(fpkm[[1]]),syms_in_pw)
cancers <- names(fpkm)
DCLs.i <- DCLs.d <- DCLres <- vector('list',length(cancers))
names(DCLs.i) <- names(DCLs.d) <- names(DCLres) <- cancers
for (i in 1:length(fpkm)) { 
  gem0 <- fpkm[[i]][kept_sym,] # Only genes of few pathways
  gem1 <- expressionBasedfilter(gem0)
  gem <- varianceBasedfilter(gem1,1e-2)
  DCGLres <- DCe(gem[,groups[[i]]==1],gem[,groups[[i]]==2],
  	link.method='percent',
  	cutoff=0.1,
  	r.method='pearson',
  	q.method='BH',
  	p=0.1,
    nbins=20,
  	figname=paste('DCGL','DCL',c('s','d'),'jpg',sep='.')
  )
  dcl <- DCLres[[i]] <- DCGLres$DCLs #subset(DCGLres$DCLs,Gene.1%in%symbol & Gene.2%in%symbol)
  #cat('All DCL:',nrow(dcl),'...')
  dcl <- subset(dcl,type=='same signed'&cor.1>0&cor.2>0)
  #cat('Same signed DCLs:',nrow(dcl),'...')
  dcl.i <- subset(dcl,cor.1<=cor.2)
  #cat('Increased DCLs:',nrow(dcl.i),'...')
  dcl.d <- subset(dcl,cor.1>cor.2)
  #cat('Decreased DCLs:',nrow(dcl.d),'...\n')
  dcl.i <- strsplit(rownames(dcl.i),',')
  dcl.d <- strsplit(rownames(dcl.d),',')
  dcl.i <- sapply(dcl.i,function(x) paste(sort(x),collapse=','))
  dcl.d <- sapply(dcl.d,function(x) paste(sort(x),collapse=','))
  DCLs.i[[i]] <- dcl.i
  DCLs.d[[i]] <- dcl.d
  DCLi <- do.call(rbind,strsplit(dcl.i,','))
  DCLd <- do.call(rbind,strsplit(dcl.d,','))
  DCLi.eg <- cbind(SYM2EG(DCLi[,1]),SYM2EG(DCLi[,2]))
  DCLd.eg <- cbind(SYM2EG(DCLd[,1]),SYM2EG(DCLd[,2]))
  DCLi.eg <- DCLi.eg[!is.na(as.numeric(DCLi.eg[,1]))&!is.na(as.numeric(DCLi.eg[,2])),]
  DCLd.eg <- DCLd.eg[!is.na(as.numeric(DCLd.eg[,1]))&!is.na(as.numeric(DCLd.eg[,2])),]
  write.table(DCLi.eg,paste0('result/',cancers[i],'.DCLi.eg.txt'),row.names=F,col.names=F,sep='\t',quote=F)
  write.table(DCLd.eg,paste0('result/',cancers[i],'.DCLd.eg.txt'),row.names=F,col.names=F,sep='\t',quote=F)
  write.table(DCLi,paste0('result/',cancers[i],'.DCLi.txt'),row.names=F,col.names=F,sep='\t',quote=F)
  write.table(DCLd,paste0('result/',cancers[i],'.DCLd.txt'),row.names=F,col.names=F,sep='\t',quote=F)
}
#save.image('dcgl.RData')
#DCLi.all <- vector('list',length(cancers))

## Merge increased/decreased DCLs across cancer types
cnt=0
for (f in dir('result',pattern='i.eg.txt',full.names=T)) {
  cnt=cnt+1
  DCLi.all[[cnt]] <- read.delim(f,head=F)
}
DCLi.all <- unique(do.call(rbind,DCLi.all))

DCLd.all <- vector('list',length(cancers))
cnt=0
for (f in dir('result',pattern='d.eg.txt',full.names=T)) {
  cnt=cnt+1
  DCLd.all[[cnt]] <- read.delim(f,head=F)
}
DCLd.all <- unique(do.call(rbind,DCLd.all))
write.table(DCLi.all,'result/DCLi.all.txt',row.names=F,col.names=F,sep='\t',quote=F)
write.table(DCLd.all,'result/DCLd.all.txt',row.names=F,col.names=F,sep='\t',quote=F)

## Create fake hsa00000 files 
pwnames <- c('FINETTI_BREAST_CANCER_KINOME_RED','MONTERO_THYROID_CANCER_POOR_SURVIVAL_UP','KUMAMOTO_RESPONSE_TO_NUTLIN_3A_DN','GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_BREAK_INDUCED_REPLICATION')
ccPWs <- matrix(NA,nr=length(pwnames),nc=2) #ConCerned PWs
for (i in 1:length(pwnames)) {
  pw <- pwnames[i]
  pwid <- paste0('hsa',10000+i)
  sym.i <- scan(paste0('data/',pw),'')
  eg.i <- SYM2EG(sym.i)
  eg.i <- as.numeric(eg.i)
  eg.i <- eg.i[!is.na(eg.i)]
  content <- data.frame(eg.i,pw,pwid)
  write.table(content,paste0('data/genes_',pwid,'.txt'),row.names=F,col.names=F,sep='\t',quote=F)
  ccPWs[i,] <- c(pwid,pw)
}
write.table(ccPWs,'ccPWs.txt',row.names=F,col.names=F,sep='\t',quote=F )
