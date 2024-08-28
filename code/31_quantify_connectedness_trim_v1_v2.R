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
  
  # Import the old way
  redeemR<-Create_redeemR(redeemR.read(path=WD,thr="S",Processed=F,rdsname = "/new.VariantsGTSummary.RDS"))
  redeemR@HomoVariants <- (redeemR@V.fitered %>% filter(totalVAF > 0.35) %>% pull(Variants)) 
  
  # New method -- trimming
  VariantsGTSummary_trim5 <- redeemR.read.trim(path=WD, "S", Processed=F,
                                               rdsname="/VariantsGTSummary.S.trim5_binom.RDS",edge_trim=5)
  
  redeemR_trim5_binom<-Create_redeemR_model(VariantsGTSummary_trim5)
  redeemR_trim5_binom<- clean_redeem(redeemR_trim5_binom,fdr = 0.05)
  redeemR_trim5_binom<-clean_redeem_removehot(redeemR_trim5_binom)


  # Append the variants that Chen manually filters out already as homoplasmic
  # so they aren't in the matrix
  redeemR@HomoVariants <- c(redeemR@HomoVariants, bad_vars_cw2)
  
  # Create matrix using redeem functions
  redeemR <- Make_matrix(redeemR)
  
  ## Filter low coverage cells and intersect with trim5
  BadCells<-subset(redeemR@CellMeta,meanCov<10)$Cell
  keep_cells <- !(rownames(redeemR@Cts.Mtx.bi) %in% BadCells)
  keep_cells2 <- intersect(rownames(redeemR_trim5_binom@Cts.Mtx.bi), rownames(redeemR@Cts.Mtx.bi))
  
  # Now compute the binary matrix under sets of circumstances
  # https://github.com/chenweng1991/redeemR/blob/master/R/BuidTree.R#L237
  Cts.Mtx.bin1 <- redeemR@Cts.Mtx.bi[keep_cells2,]
  Cts.Mtx.bin1[Cts.Mtx.bin1>=1]<-1

  # Now compute cell - cell connectivity
  # Use matrix multiplication approach from CW: https://github.com/chenweng1991/redeem_reproducibility/blob/master/Issue%231.ipynb
  cell_cell_connectivity_paper <- (Cts.Mtx.bin1 %*% t(Cts.Mtx.bin1)) >= 2
  diag(cell_cell_connectivity_paper) <- 0

  cell_cell_connectivity_any <- (Cts.Mtx.bin1 %*% t(Cts.Mtx.bin1)) >= 1
  diag(cell_cell_connectivity_any) <- 0
  
  # New method -- -2
  Cts.Mtx.bin_trim1 <- redeemR_trim5_binom@Cts.Mtx.bi[keep_cells2,!(Variants_to_underscores( colnames(redeemR_trim5_binom@Cts.Mtx.bi)) %in% redeemR@HomoVariants)]
  Cts.Mtx.bin_trim1[Cts.Mtx.bin_trim1>=1]<-1
  
  cell_cell_connectivity_trim_paper <- (Cts.Mtx.bin_trim1 %*% t(Cts.Mtx.bin_trim1)) >= 2
  diag(cell_cell_connectivity_trim_paper) <- 0
  
  cell_cell_connectivity_trim_any <- (Cts.Mtx.bin_trim1 %*% t(Cts.Mtx.bin_trim1)) >= 1
  diag(cell_cell_connectivity_trim_any) <- 0
  
  data.frame(id, 
             v1_paper = sum(cell_cell_connectivity_paper),
             v1_any = sum(cell_cell_connectivity_any),
             v2_paper = sum(cell_cell_connectivity_trim_paper),
             v2_any = sum(cell_cell_connectivity_trim_any), n_cells = length(keep_cells2)) %>%
    mutate(pct_lost_any = (1- v2_any/v1_any)*100,
           pct_lost_paper = (1- v2_paper/v1_paper)*100
    )
}

lapply(c("Young1.T1.BMMC", "Youn2.BMMC", "Old1.BMMC", "Old2.BMMC"), function(x){
  analyze_connectivity_impact(x)
}) %>% rbindlist() %>% data.frame() -> pct_df_version

saveRDS(pct_df_version, file = "../output/connectivity_lost_v1v2ReDeeM.rds")


# Now do plotting
pct_df_version <- readRDS("../output/connectivity_lost_v1v2ReDeeM.rds")

data.frame(
  id = factor(c(pct_df_version$id, pct_df_version$id), levels = c("Young1.T1.BMMC", "Youn2.BMMC", "Old1.BMMC", "Old2.BMMC")),
  version = c(rep("v1", 4), rep("v2", 4)),
  mean_c = c(pct_df_version$v1_paper/pct_df_version$n_cells,
             pct_df_version$v2_paper/pct_df_version$n_cells
             )
) %>%
  ggplot(aes(x = id, y = mean_c, fill = version)) + 
  geom_bar(stat = "identity", position = position_dodge(width=0.9), color = "black", width = 0.9) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  pretty_plot(fontsize = 7)+ L_border() +  theme(legend.position = "none") +
  labs(x = "Donor", y = "Mean connectedness per cell") -> px
cowplot::ggsave2(px, file = "../final_plots/redeem_version_connectedness.pdf", width = 3, height = 1.8)


pct_df_version <- readRDS("../output/connectivity_lost_v1v2ReDeeM.rds")

data.frame(
  id = factor(c(pct_df_version$id, pct_df_version$id), levels = c("Young1.T1.BMMC", "Youn2.BMMC", "Old1.BMMC", "Old2.BMMC")),
  version = c(rep("v1", 4), rep("v2", 4)),
  mean_c = c(pct_df_version$v1_any/pct_df_version$n_cells,
             pct_df_version$v2_any/pct_df_version$n_cells
  )
) %>%
  ggplot(aes(x = id, y = mean_c, fill = version)) + 
  geom_bar(stat = "identity", position = position_dodge(width=0.9), color = "black", width = 0.9) +
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = c("firebrick", "dodgerblue3")) +
  pretty_plot(fontsize = 7)+ L_border() +  theme(legend.position = "none") +
  labs(x = "Donor", y = "Mean # connections per cell") -> py
py
cowplot::ggsave2(py, file = "../final_plots/redeem_version_connectedness_share1var.pdf", width = 3, height = 1.8)



