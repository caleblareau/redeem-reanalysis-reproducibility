library(data.table)
library(dplyr)
library(redeemR)
library(stringr)
library(BuenColors)

# Script that calls the sensitive variants from the 
# redeem object for ease of use in downstream analyses

# update this based on your computer
# https://figshare.com/articles/dataset/ReDeeM_raw_mutation_calling/24418966/1
# Directory should have these files from the published paper
base_dir <- "/Users/lareauc/Downloads/redeem-downloaded/mito_data_redeem/"  

possible <- list.files(base_dir)

process_set <- function(id, thr_one = "S"){
  WD = paste0(base_dir, id, "/")
  redeemR<-Create_redeemR(redeemR.read(path=WD,thr=thr_one,Processed=F,rdsname = "/new.VariantsGTSummary.RDS")) # This takes several minutes
  
  redeemR.VariantsGTSummary<-CW_mgatk.read(WD,Processed =T) # This is the legacy function, the new function is redeemR.read 
  redeemR.depth<-readRDS(list.files(WD, pattern = "depth", full = TRUE)[1])  # automatically detect sample-specific name
  redeemR.Variants.feature.lst<-Vfilter_v3(InputSummary=redeemR.VariantsGTSummary,depth=redeemR.depth) 
  
  GTSummary <- redeemR@GTsummary.filtered 
  depth <- redeemR@DepthSummary
  QualifiedV <- subset(redeemR@V.fitered,HomoTag=="Hetero") 
  write.table(QualifiedV,
              file = paste0(base_dir, id, "_", thr_one, ".tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
lapply(possible, process_set)
