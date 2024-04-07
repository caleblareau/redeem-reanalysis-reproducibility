library(BuenColors)
source("../00_functions.R")

# import redeem df
redeemV_data <- fread("GSM8113080_Crispr_Mouse_Batch1.mtDNA.RawGenotypes.Total.StrandBalance.txt.gz")

# Some variants are called outside of where the fragment can be determined, making this 
# data problematic
redeemV_data %>% filter(V4 == "4287_A_T") %>% filter(grepl("5661|24231", V1))

# Nevertheless, try to identify mutations that are pervasive in different cells
n_cells <- length(unique(redeemV_data$V2))
sum_stats_df <- redeemV_data %>% 
  group_by(V2, V4) %>% summarize(obs_count = n()) %>%
  ungroup() %>% group_by(V4) %>% summarize(total_cts = sum(obs_count), pos_cells = sum(obs_count >= 0)) %>%
  arrange(desc(pos_cells)) %>% mutate(pct_cells_positive = pos_cells/n_cells)

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

