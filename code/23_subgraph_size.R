library(data.table)
library(dplyr)
library(redeemR)
library(stringr)
library(BuenColors)
library(ggtree)
library(ape)
library(phangorn)
library(ggtreeExtra)
library(igraph)
source('00_functions.R')
source("00a_modified_redeem_functions.R")

# Script that calls the sensitive variants from the 
# redeem object for ease of use in downstream analyses

compute_largest_subgraph <- function(id){
  
  WD <-  paste0(base_dir, id, ".Consensus.final/") # CL specific working directory
  
  # Import data based on reproducibility notebook
  redeemR<-Create_redeemR(redeemR.read(path=WD,thr="S",Processed=F,rdsname = "/new.VariantsGTSummary.RDS"))
  redeemR@HomoVariants <- (redeemR@V.fitered %>% filter(totalVAF > 0.35) %>% pull(Variants)) # Variants in Old 2 that are at ~40%
  
  # Append the variants that Chen manually filters out already as homoplasmic
  # so they aren't in the matrix
  redeemR@HomoVariants <- c(redeemR@HomoVariants, bad_vars_cw)
  
  ## Filter low coverage cells and compute adjacency matrix
  keep_cells<-subset(redeemR@CellMeta,meanCov>=10)$Cell
  
  # Copy object to filter edges or use only edges
  redeemR_edgeRemove<-Create_redeemR(redeemR.read_edge_partition(path=WD,thr="S",Processed=F,rdsname = "/new.VariantsGTSummary.RDS", remove_edge = TRUE))
  
  # Create matrix using redeem functions
  redeemR <- make_matrix_cl(redeemR)
  redeemR_edgeRemove <- make_matrix_cl(redeemR_edgeRemove)
  
  # create matrices
  Cts.Mtx.bin1 <- redeemR@Cts.Mtx.bi
  keep_cells <- rownames(Cts.Mtx.bin1)
  Cts.Mtx.bin2 <- redeemR_edgeRemove@Cts.Mtx.bi[keep_cells, colnames(redeemR_edgeRemove@Cts.Mtx.bi) %in% colnames(Cts.Mtx.bin1)]
  
  # export connectivity
  cell_cell_connectivity1 <- (Cts.Mtx.bin1 %*% t(Cts.Mtx.bin1)) >= 2
  diag(cell_cell_connectivity1) <- 0
  
  cell_cell_connectivity2 <- (Cts.Mtx.bin2 %*% t(Cts.Mtx.bin2)) >= 2
  diag(cell_cell_connectivity2) <- 0
  
  mc1 <- graph_from_adjacency_matrix(
    cell_cell_connectivity1
  ) %>% components(mode = "weak") 
  
  mc2 <- graph_from_adjacency_matrix(
    cell_cell_connectivity2
  ) %>% components(mode = "weak") 
  
  data.frame(
    id, 
    redeem_subgraph_pct = round(max(mc1$csize)/length(keep_cells)*100,1),
    filtered_subgraph_pct = round(max(mc2$csize)/length(keep_cells)*100,1)
  )
}

rbind(
  compute_largest_subgraph("Young1.T1.BMMC"),
  compute_largest_subgraph("Youn2.BMMC"),
  compute_largest_subgraph("Old1.BMMC"),
  compute_largest_subgraph("Old2.BMMC")
) -> subgraph_sumstats

subgraph_sumstats
