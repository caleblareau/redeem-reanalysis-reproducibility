library(data.table)
library(dplyr)
library(redeemR)
require(dplyr)
require(Matrix.utils)
source("00_functions.R")

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
  
  # Count the evidence per mutation
  mat <- redeemR@Cts.Mtx[keep_cells,]

  # Summary statistics
  
  data.frame(id, what = "redeem",
             total = sum(mat>0),
             n1 = sum(mat == 1),
             n2 = sum(mat == 2),
             n3 = sum(mat == 3))%>%
    mutate(perc1 = n1/total*100) %>%
    mutate(perc2 = n2/total*100) %>% mutate(perc3p = 100 - perc1 - perc2)
  
}

lapply(c("Young1.T1.BMMC", "Youn2.BMMC","Old1.BMMC", "Old2.BMMC"), function(x){ #, 
  analyze_connectivity_impact(x)
}) %>% rbindlist() %>% data.frame() -> pct_df1
pct_df1 <- pct_df1 

df_data <- pct_df1[,c("id", "perc1", "perc2", "perc3p")] %>%
  reshape2::melt(id.vars = c("id"))
df_data$variable <- factor(df_data$variable, levels = rev(c("perc3p", "perc2", "perc1")))
df_data$id <- factor(df_data$id, levels = (c( "Old1.BMMC", "Old2.BMMC","Young1.T1.BMMC", "Youn2.BMMC")))

pW <-  df_data %>%
  ggplot(aes(x = id, y = value, fill = variable)) + 
  geom_bar(stat = "identity", color = "black", width = 0.8, position = "dodge") + 
  scale_fill_manual(values = rev(c("lightgrey", "grey", "darkgrey"))) + 
  pretty_plot(fontsize = 7) + L_border() + 
  theme(legend.position = "none")+ 
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) + 
  labs(x = "", y = "% of cell / variant calls")
pW
cowplot::ggsave2(pW, file = "../final_plots/pct_1.pdf", width = 2.8, height = 2)


