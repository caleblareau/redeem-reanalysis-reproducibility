library(BuenColors)
source("../code/00_functions.R")

# import redeem df
redeemV_data <- fread("GSM8113080_Crispr_Mouse_Batch1.mtDNA.RawGenotypes.Total.StrandBalance.txt.gz")

# Some variants are called outside of where the fragment can be determined, making this 
# data problematic in some aspects; making a note but moving on
redeemV_data %>% filter(V4 == "4287_A_T") %>% filter(grepl("5661|24231", V1))

# Nevertheless, try to identify mutations that are pervasive in different cells
n_cells <- length(unique(redeemV_data$V2))

# Define thresholds of suspicous variants
threshold_avg <- 1.3
threshold_cells <- 20

sum_stats_df <- redeemV_data %>% 
  group_by(V2, V4) %>% summarize(obs_count = n()) %>%
  ungroup() %>% group_by(V4) %>% summarize(TotalVcount = sum(obs_count), CellN = sum(obs_count >= 0)) %>%
  arrange(desc(CellN)) %>% mutate(pct_cells_positive = CellN/n_cells) %>%
  mutate(avg_support_per_pos_cell = TotalVcount/CellN) %>%
  mutate(color = case_when( CellN > threshold_cells & avg_support_per_pos_cell < threshold_avg  ~ "redeemProblem",
                            TRUE ~ "redeemOK"))

# Make similar plot
pa <- ggplot(sum_stats_df %>% arrange(desc(color)), aes(x = CellN, y = avg_support_per_pos_cell, color = color, label = V4)) +
  geom_point(size = 0.5) + scale_x_log10() + scale_y_log10() +
  geom_hline(yintercept = threshold_avg, linetype = 2) +
  pretty_plot(fontsize = 6) + labs(x = "# cells with mutation", y = "mean # molecules / + cell") +
  scale_color_manual(values = c( "lightgrey","firebrick")) +
  theme(legend.position = "none") + L_border()

pz <- redeemV_data %>%
  filter(V4 %in% (sum_stats_df %>% filter(color == "redeemProblem") %>% pull(V4))) %>%
  make_position_df() %>%
  filter(rel_position >= 0 & rel_position <= 1) %>%
  ggplot(aes(x = rel_position)) + 
  geom_histogram(bins = 51, fill = "firebrick") + scale_y_continuous(expand = c(0,0)) + 
  labs(x = "position on molecule", y = "count of alternate allele") + pretty_plot(fontsize = 6) + L_border()
pz
get_edgeFC_mouse <- function(ggploto){
  df <-  ggplot_build(ggploto)$data[[1]]
  mean(df[c(1,2, 50,51),"y"])/mean(df[c(3:49),"y"])
}
data.frame(
  what = c("redeem"),
  value = round(c(get_edgeFC_mouse(pz)), 2)
) 

cowplot::ggsave2(cowplot::plot_grid(pa, pz, nrow = 2), 
                 width = 1.3, height = 2.3, file = paste0("../final_plots/mouse_redeem_plots.pdf"))


# Identify redeem heteroplasmic variants
top_hetero_variants_redeem <- sum_stats_df %>% filter(pct_cells_positive < 0.8) %>%
  arrange(desc(pct_cells_positive)) %>% head(60) %>% pull(V4)

# Now filter the redeemV data frame and annotate molecule position
sum_stats_redeem <- redeemV_data %>% filter(V4 %in% top_hetero_variants_redeem) %>% make_position_df()

# do KS tests on all variants
redeem_df_ks <- sum_stats_redeem %>% make_ks_test_df

# Create a set of variants with bias per caller
biased_variants_redeem <- redeem_df_ks %>% filter(statistic > 0.35) %>% pull(variant)

# Function to plot positions of alt alleles on molecule
make_plot <- function(var_df, var_levels, color, title) {
  n_vars = length(unique(var_df$variant))
  pP <- ggplot(var_df %>%  
                 filter(rel_position >= 0 & rel_position <= 1) %>%
                 mutate(variant2 = factor(as.character(variant), levels = var_levels)), 
               aes( x = rel_position)) +
    geom_histogram(color = color, fill = color, bins = 30) + 
    facet_wrap(~variant2, scales = "free_y") +
    ggtitle(paste0("Mouse_Batch1", " - ",title," top variants")) + theme_bw() +
    labs(x = "Position on molecule", y = "abundance of alt allele") 
  
  ggsave(pP, file = paste0("Mouse_Batch1_",title,".png"),
         width = (n_vars)^(1/3)*4, height = (n_vars)^(1/3)*2, dpi = 200)
  title
}

pRlow <- make_plot(
  sum_stats_redeem %>% filter( !(variant %in% biased_variants_redeem)),
  top_hetero_variants_redeem, "dodgerblue3", "Redeem_KSretained")

pRhigh <- make_plot(
  sum_stats_redeem %>% filter(variant %in% biased_variants_redeem),
  top_hetero_variants_redeem, "firebrick", "Redeem_KSremoved")

