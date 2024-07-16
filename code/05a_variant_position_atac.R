library(BuenColors)
source("00_functions.R")

sample_id <- "Young2.BMMC"
redeem_id <- "Youn2.BMMC"

redeemV_data <- fread(paste0("../data/arc_data_from_atac/",sample_id,".ATAC.RawGenotypes.Sensitive.gz"))

# Load scatters
redeem_df <- fread(paste0("../data/redeem-sensitive-calls/",redeem_id,".Consensus.final_S.tsv"))

# Define thresholds of suspicous variants
threshold_avg <- 1.3
threshold_cells <- 20

sens_df_annotated <- redeem_df %>%
  filter(!c(Variants %in% bad_vars_cw)) %>%
  mutate(avg_support_per_pos_cell = TotalVcount/CellN) %>%
  mutate(color = case_when(CellN > threshold_cells & avg_support_per_pos_cell < threshold_avg  ~ "redeemProblem",
                           TRUE ~ "redeemOK"))

p2 <- redeemV_data %>%
  filter(V4 %in% (sens_df_annotated %>% filter(color == "redeemProblem") %>% pull(Variants))) %>%
  make_position_df() %>%
  ggplot(aes(x = rel_position)) + 
  geom_histogram(bins = 51, fill = "red") + scale_y_continuous(expand = c(0,0)) + 
  labs(x = "position on molecule", y = "count of alternate allele") + pretty_plot(fontsize = 6) + L_border()
p2

cowplot::ggsave2(p2, 
                 width = 1.4, height = 1.4, file = paste0("../final_plots/", sample_id, "_atac_plot.pdf"))


get_edgeFC <- function(ggploto){
  df <-  ggplot_build(ggploto)$data[[1]]
  mean(df[c(1:3, 49,50,51),"y"])/mean(df[c(4:48),"y"])
}
data.frame(
  what = c( "redeem_low"),
  value = round(c(get_edgeFC(p2)), 2)
) 

