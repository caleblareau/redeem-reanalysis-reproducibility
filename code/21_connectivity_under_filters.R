library(data.table)
library(dplyr)
library(redeemR)
library(stringr)
library(BuenColors)
library(ggtree)
library(ape)
library(phangorn)
library(ggtreeExtra)
source('00_functions.R')
source("00a_modified_redeem_functions.R")


connectivity_recompute <- function(id){
  
  WD <-  paste0(base_dir, id, ".Consensus.final/") # CL specific working directory
  
  # Import data based on reproducibility notebook
  redeemR<-Create_redeemR(redeemR.read(path=WD,thr="S",Processed=F,rdsname = "/new.VariantsGTSummary.RDS"))
  redeemR@HomoVariants <- (redeemR@V.fitered %>% filter(totalVAF > 0.35) %>% pull(Variants)) # Variants in Old 2 that are at ~40%
  
  # Append the variants that Chen manually filters out already as homoplasmic
  # so they aren't in the matrix
  redeemR@HomoVariants <- c(redeemR@HomoVariants, bad_vars_cw)
  
  ## Filter low coverage cells and compute adjacency matrix
  keep_cells<-subset(redeemR@CellMeta,meanCov>=10)$Cell
  
  # Copy object to filter edges with two different ways
  redeemR_edgeRemove<-Create_redeemR(redeemR.read_edge_partition(path=WD,thr="S",Processed=F,rdsname = "/new.VariantsGTSummary.RDS", remove_edge = TRUE, ratio_based_filter = TRUE))
  redeemR_edgeRemove_cw<-Create_redeemR(redeemR.read_edge_partition(path=WD,thr="S",Processed=F,rdsname = "/new.VariantsGTSummary.RDS", remove_edge = TRUE, ratio_based_filter = FALSE))
  
  # Create matrix using redeem functions
  redeemR <- make_matrix_cl(redeemR)
  redeemR_edgeRemove <- make_matrix_cl(redeemR_edgeRemove)
  redeemR_edgeRemove_cw <- make_matrix_cl(redeemR_edgeRemove_cw)
  
  # create matrices
  Cts.Mtx.bin1 <- redeemR@Cts.Mtx.bi
  keep_cells <- rownames(Cts.Mtx.bin1)
  Cts.Mtx.bin2 <- redeemR_edgeRemove@Cts.Mtx.bi[keep_cells, colnames(redeemR_edgeRemove@Cts.Mtx.bi) %in% colnames(Cts.Mtx.bin1)]
  Cts.Mtx.bin3 <- redeemR_edgeRemove_cw@Cts.Mtx.bi[keep_cells, colnames(redeemR_edgeRemove_cw@Cts.Mtx.bi) %in% colnames(Cts.Mtx.bin1)]
  
  # export connectivity
  cell_cell_connectivity1 <- (Cts.Mtx.bin1 %*% t(Cts.Mtx.bin1)) >= 2
  diag(cell_cell_connectivity1) <- 0
  
  cell_cell_connectivity2 <- (Cts.Mtx.bin2 %*% t(Cts.Mtx.bin2)) >= 2
  diag(cell_cell_connectivity2) <- 0
  
  cell_cell_connectivity3 <- (Cts.Mtx.bin3 %*% t(Cts.Mtx.bin3)) >= 2
  diag(cell_cell_connectivity3) <- 0
  
  sum_stats_df <-   data.frame(id, 
                               connections_original = sum(cell_cell_connectivity1),
                               connections_edgeRemove = sum(cell_cell_connectivity2),
                               connections_edgeRemove_cw = sum(cell_cell_connectivity3)
  )
  sum_stats_df
}

stat_df <- rbind(
  connectivity_recompute("Young1.T1.BMMC"),
  connectivity_recompute("Youn2.BMMC"),
  connectivity_recompute("Old2.BMMC"),
  connectivity_recompute("Old1.BMMC")
) %>%
  mutate(perc_lost_CL = (1-connections_edgeRemove/connections_original) * 100)%>%
  mutate(perc_lost_CW = (1-connections_edgeRemove_cw/connections_original) * 100)
stat_df

