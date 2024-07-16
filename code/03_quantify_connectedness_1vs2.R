library(data.table)
library(dplyr)
library(redeemR)
require(dplyr)
require(Matrix.utils)
source("00_functions.R")

# example
id = "Young1.T1.BMMC"

id <- "Youn2.BMMC"
analyze_connectivity_impact <- function(id){
  print(id)
  WD <-  paste0(base_dir, id, ".Consensus.final/") # CL specific working directory
  
  # Import data based on reproducibility notebook
  redeemR<-Create_redeemR(redeemR.read(path=WD,thr="S",Processed=F,rdsname = "/new.VariantsGTSummary.RDS"))
  redeemR@HomoVariants <- (redeemR@V.fitered %>% filter(totalVAF > 0.35) %>% pull(Variants)) # Variants in Old 2 that are at 40%
  
  (redeemR@V.fitered %>% filter(totalVAF > 0.35) %>% filter(HomoTag != "Homo"))
  
  # Append the variants that Chen manually filters out already as homoplasmic
  # so they aren't in the matrix
  redeemR@HomoVariants <- c(redeemR@HomoVariants, bad_vars_cw)
  
  # Create matrix using redeem functions
  redeemR <- Make_matrix(redeemR)
  
  
  ## Filter low coverage cells
  BadCells<-subset(redeemR@CellMeta,meanCov<10)$Cell
  keep_cells <- !(rownames(redeemR@Cts.Mtx.bi) %in% BadCells)

  # Now compute the binary matrix under two sets of circumstances
  # https://github.com/chenweng1991/redeemR/blob/master/R/BuidTree.R#L237

  Cts.Mtx.bin1 <- redeemR@Cts.Mtx[keep_cells,]
  Cts.Mtx.bin1[Cts.Mtx.bin1>=1]<-1

  Cts.Mtx.bin2 <- redeemR@Cts.Mtx[keep_cells,]
  Cts.Mtx.bin2[Cts.Mtx.bin2 <= 1]<-0
  Cts.Mtx.bin2[Cts.Mtx.bin2 >= 2]<-1
  
  # Now compute cell - cell connectivity
  # Use matrix multiplication approach from CW: https://github.com/chenweng1991/redeem_reproducibility/blob/master/Issue%231.ipynb

  cell_cell_connectivity1 <- (Cts.Mtx.bin1 %*% t(Cts.Mtx.bin1)) >= 2
  diag(cell_cell_connectivity1) <- 0
  
  cell_cell_connectivity2 <- (Cts.Mtx.bin2 %*% t(Cts.Mtx.bin2)) >= 2
  diag(cell_cell_connectivity2) <- 0

  cell_cell_connectivity3 <- (Cts.Mtx.bin1 %*% t(Cts.Mtx.bin1)) >= 1
  diag(cell_cell_connectivity3) <- 0
  
  # Summarize what happens
  pct_connections_lost <- (sum(cell_cell_connectivity1) - sum(cell_cell_connectivity2))/sum(cell_cell_connectivity1)*100
  
  data.frame(id, 
             original_connections = sum(cell_cell_connectivity1),
             connections_after_filter1 = sum(cell_cell_connectivity2),
             connections_wOnly1 = sum(cell_cell_connectivity3), n_cells = sum(keep_cells),
             pct_connections_lost = round(pct_connections_lost, 2))
}

lapply(c("Young1.T1.BMMC", "Youn2.BMMC", "Old1.BMMC", "Old2.BMMC"), function(x){
  analyze_connectivity_impact(x)
}) %>% rbindlist() %>% data.frame() -> pct_df12

(pct_df12 %>%
  mutate(original = original_connections/n_cells, 
         only1 = connections_wOnly1/n_cells))[,c("id", "original", "only1")]


(pct_df12 %>%
    mutate(original = original_connections/n_cells/n_cells*100, 
           only1 = connections_wOnly1/n_cells/n_cells*100))[,c("id", "original", "only1")]


pX <- pct_df12 %>%
  ggplot(aes(x = id, y = pct_connections_lost)) + 
  geom_bar(stat = "identity", color = "black", fill = "lightgrey", width = 0.6) + 
  pretty_plot(fontsize = 7) + L_border() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) + 
  labs(x = "", y = "% connectivity lost")
cowplot::ggsave2(pX, file = "../final_plots/connectivity_lost.pdf", width = 1.4, height = 1)


pct_go <-  pct_df12 %>% mutate(Weng = original_connections/n_cells, 
                               zPostFilter = connections_after_filter1/n_cells)

pY <- pct_go[,c("id","Weng", "zPostFilter")] %>%
  reshape2::melt(id.vars = "id") %>% 
  ggplot(aes(x = id, y = value, fill = variable)) + 
  geom_bar(stat = "identity", color = "black", width = 0.6, position = position_dodge2()) + 
  pretty_plot(fontsize = 7) + L_border() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 70)) + 
  labs(x = "", y = "mean connectivity") + 
  scale_fill_manual(values = c("grey1", "grey")) +
  theme(legend.position = "none")
cowplot::ggsave2(pY, file = "../final_plots/n_connections_more1.pdf", width = 2, height = 1)

#############

boosted_go <-  pct_df12 %>% mutate(Weng = original_connections/n_cells, 
                               zPostMod = connections_wOnly1/n_cells)

pZ <- boosted_go[,c("id","Weng", "zPostMod")] %>%
  reshape2::melt(id.vars = "id") %>% 
  ggplot(aes(x = id, y = value, fill = variable)) + 
  geom_bar(stat = "identity", color = "black", width = 0.6, position = position_dodge2()) + 
  pretty_plot(fontsize = 7) + L_border() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 400)) + 
  labs(x = "", y = "mean connectivity") + 
  scale_fill_manual(values = c("grey1", "grey")) +
  theme(legend.position = "none")
pZ
cowplot::ggsave2(pZ, file = "../final_plots/n_connections_only1.pdf", width = 2, height = 1)




