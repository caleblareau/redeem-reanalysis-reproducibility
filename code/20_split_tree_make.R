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

# Script that calls the sensitive variants from the 
# redeem object for ease of use in downstream analyses

partition_tree_edges <- function(id){
  
  WD <-  paste0(base_dir, id, ".Consensus.final/") # CL specific working directory
  
  # Import data based on reproducibility notebook
  redeemR<-Create_redeemR(redeemR.read(path=WD,thr="S",Processed=F,rdsname = "/new.VariantsGTSummary.RDS"))
  redeemR@HomoVariants <- (redeemR@V.fitered %>% filter(totalVAF > 0.35) %>% pull(Variants)) # Variants in Old 2 that are at ~40%
  
  # Append the variants that Chen manually filters out already as homoplasmic
  # so they aren't in the matrix
  redeemR@HomoVariants <- c(redeemR@HomoVariants, bad_vars_cw)
  
  ## Filter low coverage cells and compute adjacency matrix
  keep_cells<-subset(redeemR@CellMeta,meanCov>=10)$Cell
  
  # Create matrix using redeem functions
  redeemR <- Make_matrix(redeemR)
  bad_vars_cw %in% colSums(redeemR@Cts.Mtx.bi) 
  
  ## Filter low coverage cells and compute adjacency matrix
  BadCells<-subset(redeemR@CellMeta,meanCov<10)$Cell
  keep_cells <- !(rownames(redeemR@Cts.Mtx.bi) %in% BadCells)
  
  # now filter out the instances where 
  # we only have 1 molecule supporting the interaction
  Cts.Mtx.bin1 <- redeemR@Cts.Mtx[keep_cells,]
  Cts.Mtx.bin1[Cts.Mtx.bin1>=1]<-1
  
  # number 2 is only 2+ reads
  Cts.Mtx.bin2 <- redeemR@Cts.Mtx[keep_cells,]
  Cts.Mtx.bin2[Cts.Mtx.bin2 <= 1]<-0
  Cts.Mtx.bin2[Cts.Mtx.bin2 >= 2]<-1
  
  # number 3 is only == 1
  Cts.Mtx.bin3 <- redeemR@Cts.Mtx[keep_cells,]
  Cts.Mtx.bin3[Cts.Mtx.bin3 >= 2]<-0
  
  # export connectivity
  cell_cell_connectivity1 <- (Cts.Mtx.bin1 %*% t(Cts.Mtx.bin1)) >= 2
  diag(cell_cell_connectivity1) <- 0
  
  cell_cell_connectivity2 <- (Cts.Mtx.bin2 %*% t(Cts.Mtx.bin2)) >= 2
  diag(cell_cell_connectivity2) <- 0
  
  cell_cell_connectivity3 <- (Cts.Mtx.bin3 %*% t(Cts.Mtx.bin3)) >= 2
  diag(cell_cell_connectivity3) <- 0
 

  # create a data frame sumamarizing connectivity
  sum_stats_df <-   data.frame(id, 
                               connections_original = sum(cell_cell_connectivity1),
                               connections_2p = sum(cell_cell_connectivity2),
                               connections_only1 = sum(cell_cell_connectivity3)
  )
  saveRDS(sum_stats_df, file = paste0("../output/trees_partition/connectivity_all_1_2p_", id, "_df.rds"))
  
  
  # function to make the weights for their jaccard thingy
  make_weight_df <- function(Cts.bin){
    # Pull out weights fromt the data object; update with N/As with 1 as described in manuscript
    weightdf<-data.frame(Variants=colnames(Cts.bin)) %>% merge(.,V.weight,by="Variants",all.x = T,sort = F) 
    weight <- weightdf %>% pull(weight)
    weight[is.na(weight)] <- 1
    weight
  }
  
  # do the weights
  weight1 <- make_weight_df(Cts.Mtx.bin1)
  weight2 <- make_weight_df(Cts.Mtx.bin2)
  weight3 <- make_weight_df(Cts.Mtx.bin3)

  # Compute weighted jaccard distance metric in the paper
  w_j1 <- quick_w_jaccard_cl(Cts.Mtx.bin1,weight1) %>% data.matrix()
  w_j2 <- quick_w_jaccard_cl(Cts.Mtx.bin2,weight2) %>% data.matrix()
  w_j3 <- quick_w_jaccard_cl(Cts.Mtx.bin3,weight3) %>% data.matrix()

  # Make phylogenetic tree
  phylo1 <- nj(w_j1)
  phylo2 <- nj(w_j2)
  phylo3 <- nj(w_j3)

  save(phylo1, phylo2, phylo3,  file = paste0("../output/trees_partition/trees_all_1_2p_", id, ".rda"))
  
  id
}

partition_tree_edges("Young1.T1.BMMC")
partition_tree_edges("Youn2.BMMC")
partition_tree_edges("Old2.BMMC")
partition_tree_edges("Old1.BMMC")

readRDS(paste0("../output/trees_partition/connectivity_all_1_2p_Young1.T1.BMMC_df.rds"))
