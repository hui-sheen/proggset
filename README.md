# A Novel Strategy to Identify Prognosis Relevant Gene Sets in Cancers
## Step 1: Single-gene survival analysis
#Check code/Step1.survival.R  
## Step 2: GSEA based on single gene survival analysis results
#Check code/Step2.fgsea.R  
## Step 3: Gene set survival analysis based on CGES
#Check code/Step3.gsea_flow2_CGES.R
## Step 4: Crosstalk analysis using DCGL & CSPN
#Check code/Step4.crosstalk.R
#Unzip cspn  
cd cspn  
#Deploy three aspects of files in cspn  
##background gene-gene network @ cspn/reference/PPI_network/  
##Gene set definition identified with EntrezGeneIDs @ cspn/reference/original_gene_sets/  
##Concerned pathway list @ cspn/reference/concerned_pathway_list/  
#$nw is project keyword  
#$nwFile designates file name for background network (say, DCLi.all.txt or DCLd.all.txt)  
perl ./cspn.pl $nwFile ./reference/concerned_pathway_list/sigPWs.27.lst ./result/$nw.ctalk ./result/$nw.clink  
