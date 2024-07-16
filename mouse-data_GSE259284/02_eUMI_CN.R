library(BuenColors)
library(data.table)
source("../code/00_functions.R")

# import redeem df
redeemV_data <- fread("GSM8113080_Crispr_Mouse_Batch1.mtDNA.RawGenotypes.Total.StrandBalance.txt.gz")
dt <- fread(paste0("../../redeem-downloaded/mito_data_redeem/Young1.T1.BMMC.Consensus.final/RawGenotypes.Sensitive.StrandBalance"))
summary(redeemV_data$V8)
summary(dt$V8)

data.frame(
  what = c(rep("mouse", dim(redeemV_data)[1]),
           rep("human", dim(dt)[1])),
  eUMI_count = c(redeemV_data[["V8"]],dt[["V8"]])
) -> long_df
long_df%>%
  ggplot(aes(x = eUMI_count, color = what )) + 
  stat_ecdf() + 
  scale_x_log10() + 
  coord_cartesian(xlim = c(1, 30)) +
  labs(x = "# PCR duplicates / eUMI", y = "Cumulative frequency") +
  pretty_plot(fontsize = 7) + L_border() + theme(legend.position = "none") -> px
cowplot::ggsave2(px, file = "../final_plots/mouse_human_ecdf.pdf", width = 1.5, height = 1.5)
table(long_df$what)
mean(redeemV_data[["V8"]] == 1)
