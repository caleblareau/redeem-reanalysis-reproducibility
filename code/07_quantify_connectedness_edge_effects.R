library(data.table)
library(dplyr)
library(redeemR)

source("00_functions.R")

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
  bad_vars_cw <- c("310_T_C","3109_T_C","309_C_T", "5764_C_T") # variants that we agree are bad / already filtered; pretend these are homoplasmic to remove
  redeemR@HomoVariants <- c(redeemR@HomoVariants, bad_vars_cw)
  
  # Identify potential problematic variants from the KS test
  redeemV_data <- fread(paste0(WD, "/RawGenotypes.Sensitive.StrandBalance"))
  redeem_called_vars <- redeemR@V.fitered$Variants
  sum_stats_redeem <- redeemV_data %>% filter(V4 %in% redeem_called_vars) %>% 
    make_position_df() %>% make_ks_test_df()
  suspected_bad <- sum_stats_redeem %>% filter(statistic > 0.35) %>% pull(variant)
  suspected_bad_formatted <- paste0("Variants", gsub("_", "", suspected_bad))
  
  # Create matrix using redeem functions
  redeemR <- Make_matrix(redeemR)
  bad_vars_cw %in% colSums(redeemR@Cts.Mtx.bi) 
  
  ## Filter low coverage cells
  BadCells<-subset(redeemR@CellMeta,meanCov<10)$Cell
  keep_cells <- !(rownames(redeemR@Cts.Mtx.bi) %in% BadCells)
  Cts.Mtx.bin <- redeemR@Cts.Mtx.bi[keep_cells,]
  
  # Now compute cell - cell connectivity
  # 0 out self connectivity
  # Use matrix multiplication approach from CW: https://github.com/chenweng1991/redeem_reproducibility/blob/master/Issue%231.ipynb
  cell_cell_connectivity <- (Cts.Mtx.bin %*% t(Cts.Mtx.bin)) >= 2
  diag(cell_cell_connectivity) <- 0
  
  # now filter out the position biased variants
  boo_keep <- !(colnames(Cts.Mtx.bin) %in% suspected_bad_formatted)
  sum(boo_keep)
  
  # recompute the connectivity
  # Using matrix multiplication approach from CW: https://github.com/chenweng1991/redeem_reproducibility/blob/master/Issue%231.ipynb
  Cts.Mtx.bin.filtered <- Cts.Mtx.bin[,boo_keep ]
  cell_cell_connectivity_post_filter <- (Cts.Mtx.bin.filtered %*% t(Cts.Mtx.bin.filtered)) >= 2
  diag(cell_cell_connectivity_post_filter) <- 0
  
  # Summarize what happens
  pct_connections_lost <- (sum(cell_cell_connectivity) - sum(cell_cell_connectivity_post_filter))/sum(cell_cell_connectivity)*100
  
  data.frame(id, 
             original_connections = sum(cell_cell_connectivity),
             connections_after_filter = sum(cell_cell_connectivity_post_filter),
             pct_connections_lost = round(pct_connections_lost, 2))
}

lapply(c("Young1.T1.BMMC", "Youn2.BMMC", "Old1.BMMC", "Old2.BMMC"), function(x){
  analyze_connectivity_impact(x)
}) %>% rbindlist() %>% data.frame() -> pct_df_ks_lost

pct_df_ks_lost



