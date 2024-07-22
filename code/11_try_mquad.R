library(data.table)
library(dplyr)
library(redeemR)
require(dplyr)
require(Matrix.utils)
library(Matrix)
source("00_functions.R")

# Modified function to get the depth per variant / cell
Add_DepthMatrix_CL <- function(object,QualifiedTotalCtsFile=""){
  qtc <- fread(QualifiedTotalCtsFile)
  colnames(qtc)<-c("Cell","Pos","T","LS","S","VS")
  Dic<-gsub("Variants","",colnames(object@Cts.Mtx.bi)) %>% substr(.,1,nchar(.)-2) %>% as.integer %>% data.frame(Variants=colnames(object@Cts.Mtx.bi),Pos=.)
  QualifiedTotalCts.subset<-qtc[qtc$Cell %in% row.names(object@Cts.Mtx.bi),] %>% merge(.,Dic,by="Pos") 
  DepthMatrix<-dcast.data.table(QualifiedTotalCts.subset[,c("Cell", "Variants","S")],Cell~Variants) %>% data.frame() %>% tibble::column_to_rownames("Cell") %>% as.matrix
  if (all(dim(object@Cts.Mtx.bi)==dim(DepthMatrix))){
    object@Ctx.Mtx.depth<-DepthMatrix[row.names(object@Cts.Mtx.bi),colnames(object@Cts.Mtx.bi)]
  }else{
    print(dim(object@Cts.Mtx.bi))
    print(dim(DepthMatrix))
    print("Check the input QualifiedTotalCts, the dimension doesn't match")
  }
  return(object)
}

# example
id = "Young1.T1.BMMC"

analyze_connectivity_impact <- function(id){
  print(id)
  WD <-  paste0(base_dir, id, ".Consensus.final/") # CL specific working directory
  
  # Import data based on reproducibility notebook
  redeemR<-Create_redeemR(redeemR.read(path=WD,thr="S",Processed=F,rdsname = "/new.VariantsGTSummary.RDS"))
  redeemR@HomoVariants <- (redeemR@V.fitered %>% filter(totalVAF > 0.35) %>% pull(Variants)) # Variants in Old 2 that are at 40%
  
  # Append the variants that Chen manually filters out already as homoplasmic
  # so they aren't in the matrix
  redeemR@HomoVariants <- c(redeemR@HomoVariants, bad_vars_cw)
  
  # Create matrix using redeem functions
  redeemR <- Make_matrix(redeemR)
  
  ## Filter low coverage cells
  BadCells<-subset(redeemR@CellMeta,meanCov<10)$Cell
  keep_cells <- !(rownames(redeemR@Cts.Mtx.bi) %in% BadCells)
  redeemR <- Add_DepthMatrix_CL(redeemR, QualifiedTotalCtsFile = paste0(WD, "QualifiedTotalCts"))
  
  writeMM(t(redeemR@Cts.Mtx), file = paste0("../../redeem-downloaded/large_for_mquad/", id, ".AD.mtx"))
  writeMM(t(Matrix(redeemR@Ctx.Mtx.depth, sparse = TRUE)), file = paste0("../../redeem-downloaded/large_for_mquad/", id, ".DP.mtx"))
  write.table(data.frame(colnames(redeemR@Ctx.Mtx.depth)), file = paste0("../../redeem-downloaded/large_for_mquad/", id, ".variants.out.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

lapply(c( "Youn2.BMMC", "Old1.BMMC", "Old2.BMMC"), function(x){ # "Young1.T1.BMMC",
  analyze_connectivity_impact(x)
}) %>% rbindlist() %>% data.frame() -> pct_df12

# Next, run MQuad with default parameters to assess what gets called variant
#/Users/lareauc/Library/Python/3.9/bin/mquad -m Young1.T1.BMMC.AD.mtx,Young1.T1.BMMC.DP.mtx -o Young1.T1.BMMC -p 16
#/Users/lareauc/Library/Python/3.9/bin/mquad -m Youn2.BMMC.AD.mtx,Youn2.BMMC.DP.mtx -o Youn2.BMMC -p 16
#/Users/lareauc/Library/Python/3.9/bin/mquad -m Old1.BMMC.AD.mtx,Old1.BMMC.DP.mtx -o Old1.BMMC -p 16
#/Users/lareauc/Library/Python/3.9/bin/mquad -m Old2.BMMC.AD.mtx,Old2.BMMC.DP.mtx -o Old2.BMMC -p 16

