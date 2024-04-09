library(data.table)
library(dplyr)
library(redeemR)
library(matrixStats)
library(ggbeeswarm)
source("00_functions.R")

# Function modified from here: https://github.com/chenweng1991/redeemR/blob/93332c648b5cd9310609271da3a199e9e9f98167/R/BuidTree.R#L919
quick_w_jaccard_cl<-function(M,w){ 
  total<-M %*% w
  a<-M %*% (Matrix::t(M)*w)
  b<-as.numeric(total) - a
  c<-Matrix::t(b)
  disimilarity<-round(1-a/(a+b+c+0.00001),4)
  #distance<-as.dist(disimilarity)
  
  # One other modification to force same cell to be its closest neighbor always
  diag(disimilarity) <- -0.01
  return(disimilarity)
}

# same function as here for comparison: 
# https://github.com/chenweng1991/redeemR/blob/93332c648b5cd9310609271da3a199e9e9f98167/R/BuidTree.R#L919
quick_w_jaccard<-function(M,w){ 
  total<-M %*% w
  a<-M %*% (Matrix::t(M)*w)
  b<-as.numeric(total) - a
  c<-Matrix::t(b)
  disimilarity<-round(1-a/(a+b+c+0.00001),4)
  #distance<-as.dist(disimilarity)
  return(disimilarity)
}

analyze_neighbors_impact <- function(id){
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
  
  ## Filter low coverage cells and compute adjacency matrix
  BadCells<-subset(redeemR@CellMeta,meanCov<10)$Cell
  keep_cells <- !(rownames(redeemR@Cts.Mtx.bi) %in% BadCells)
  Cts.Mtx.bin <- redeemR@Cts.Mtx.bi[keep_cells,]
  
  # now filter out the position biased variants
  boo_keep <- !(colnames(Cts.Mtx.bin) %in% suspected_bad_formatted)
  sum(boo_keep)
  Cts.Mtx.bin.filtered <- Cts.Mtx.bin[,boo_keep ]
  
  # Pull out weights fromt the data object; update with N/As with 1 as described in manuscript
  weightdf<-data.frame(Variants=colnames(Cts.Mtx.bin)) %>% merge(.,V.weight,by="Variants",all.x = T,sort = F) 
  weight <- weightdf %>% pull(weight)
  weight[is.na(weight)] <- 1
  
  # Compute weighted jaccard distance metric in the paper
  w_j <- quick_w_jaccard_cl(Cts.Mtx.bin,weight)
  w_j_cw <- quick_w_jaccard(Cts.Mtx.bin,weight)
  
  # one difference is a small pseudocount for divisions by zero
  stopifnot(sum(w_j < 1) ==  sum(w_j_cw < 1)) 
  
  # Define the rank per cell with rowRanks, breaking ties with random placement
  rr_original <- rowRanks(data.matrix(w_j), ties.method = "random")
  
  # Other difference is to make sure that the same cell is always its own neighbor
  # about 100 instances where a cell is tied with its neighbor and randomly placed otherwise
  stopifnot(all(diag(rr_original) ==1))
  
  w_j2 <- quick_w_jaccard_cl(Cts.Mtx.bin.filtered,weight[boo_keep])
  rr_updated <- rowRanks(data.matrix(w_j2), ties.method = "random")
  
  # Now compute
  compute_same_nn_pre_post_variant_filtering <- function(n_neighbors){
    same_pre_post_per_cell <- rowSums((rr_updated <= (n_neighbors + 1)) * 
                                        (rr_original <= (n_neighbors + 1))) - 1
    round((1 - mean(same_pre_post_per_cell/n_neighbors))*100, 2)
  }
  
  data.frame(id,
             nn_01 = compute_same_nn_pre_post_variant_filtering(1),
             nn_05 = compute_same_nn_pre_post_variant_filtering(5),
             nn_10 = compute_same_nn_pre_post_variant_filtering(10),
             nn_50 = compute_same_nn_pre_post_variant_filtering(50)
  )
}

lapply(c("Young1.T1.BMMC", "Youn2.BMMC", "Old1.BMMC", "Old2.BMMC"), function(x){
  analyze_neighbors_impact(x)
}) %>% rbindlist() %>% data.frame() -> pct_nn_lost

p1 <- ggplot(pct_nn_lost %>% reshape2::melt(id.vars = "id"), aes(x = variable, y = value, color = id)) +
  geom_point() +
  labs(x = "", y = "% nearest neighbors changed") + 
  pretty_plot(fontsize = 8) + L_border() + 
  scale_y_continuous(limits = c(0,100)) +
  geom_hline(yintercept = c(0, 100), linetype = 2) + 
  scale_color_manual(values = jdb_palette("corona")[4:1]) +
  theme(legend.position = "none")

cowplot::ggsave2(p1, file = "../final_plots/nn_results.pdf", width = 1.9, height = 1.7)
