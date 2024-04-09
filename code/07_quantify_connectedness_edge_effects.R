library(data.table)
library(dplyr)
library(redeemR)
library(matrixStats)

source("00_functions.R")

# example
id = "Young1.T1.BMMC"

analyze_connectivity_impact <- function(id){
  print(id)
  WD <-  paste0(base_dir, id, ".Consensus.final/") # CL specific working directory
  
  # Import data based on reproducibility notebook
  redeemR<-Create_redeemR(redeemR.read(path=WD,thr="S",Processed=F,rdsname = "/new.VariantsGTSummary.RDS"))
  redeemR@HomoVariants <- (redeemR@V.fitered %>% filter(totalVAF > 0.35) %>% pull(Variants)) # Variants in Old 2 that are at ~40%
  
  # Append the variants that Chen manually filters out already as homoplasmic
  # so they aren't in the matrix
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
             pct_connections_lost = round(pct_connections_lost, 2),
             median_connections_pre = median(colSums(data.matrix(cell_cell_connectivity))) %>% round(1),
             median_connections_post = median(colSums(data.matrix(cell_cell_connectivity_post_filter))) %>% round(1),
             mean_connections_pre = mean(colSums(data.matrix(cell_cell_connectivity))) %>% round(1),
             mean_connections_post = mean(colSums(data.matrix(cell_cell_connectivity_post_filter))) %>% round(1)
             )
}

lapply(c("Young1.T1.BMMC", "Youn2.BMMC", "Old1.BMMC", "Old2.BMMC"), function(x){
  analyze_connectivity_impact(x)
}) %>% rbindlist() %>% data.frame() -> pct_df_ks_lost

pct_df_ks_lost

pX <- pct_df_ks_lost %>%
  ggplot(aes(x = id, y = pct_connections_lost)) + 
  geom_bar(stat = "identity", color = "black", fill = "lightgrey", width = 0.6) + 
  pretty_plot(fontsize = 7) + L_border() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) + 
  labs(x = "", y = "% connectivity lost")
cowplot::ggsave2(pX, file = "../final_plots/connectivity_lost_ks.pdf", width = 1.7, height = 1.7)

pct_df_ks_lost[,c("id", "mean_connections_pre","mean_connections_post")] %>%
  reshape2::melt(id.vars = "id") %>%
  ggplot(aes(x = id, fill = variable, y = value)) + 
  geom_bar(stat = "identity", position=position_dodge()) +
  scale_fill_manual(values = c("lightgrey", "darkgrey")) +
  pretty_plot(fontsize = 8)  + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position = "none") -> p1
cowplot::ggsave2(p1, file = "../final_plots/pre_post_mean.pdf", width = 3.5, height = 1.7)
