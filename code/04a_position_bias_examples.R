library(BuenColors)
source("00_functions.R")

# Main Function that processes summary counts
process_one_sample <- function(sample_id){
  
  # import redeem df
  redeemV_data <- fread(paste0(base_dir,sample_id,".Consensus.final/RawGenotypes.Sensitive.StrandBalance"))
  
  # import and select top 60 variants per variant caller
  sens_call_df <- fread(paste0("../data/redeem-sensitive-calls/",sample_id,".Consensus.final_S.tsv"))
  top_hetero_variants_redeem <- sens_call_df %>% filter(totalVAF < 0.8) %>%
    arrange(desc(CellN)) %>% head(60) %>% pull(Variants)
  
  ## mgatk calls 
  mgatk_df_in <- fread(paste0("../data/mgatk-varstats/D_",sample_id,"_hg38_v20-mtMask.variant_stats.tsv.gz"))
  mgatk_df <- mgatk_df_in %>% filter(strand_correlation > 0.65 & n_cells_conf_detected >= 3 & (log10(vmr) > -2) & mean < 0.8)
  top_hetero_variants_mgatk <- mgatk_df %>% arrange(desc(n_cells_conf_detected)) %>% 
    mutate(redeem_formatted = paste0(position, "_", gsub(">", "_", nucleotide))) %>%
    pull(redeem_formatted) %>% head(60)
  
  # Now filter the redeemV data frame and annotate molecule position
  sum_stats_redeem <- redeemV_data %>% filter(V4 %in% top_hetero_variants_redeem) %>% make_position_df()
  sum_stats_mgatk <- redeemV_data %>% filter(V4 %in% top_hetero_variants_mgatk) %>% make_position_df()
  
  # do KS tests on all variants
  redeem_df_ks <- sum_stats_redeem %>% make_ks_test_df
  mgatk_df_ks <- sum_stats_mgatk %>% make_ks_test_df
  
  # Create a set of variants with bias per caller
  biased_variants_redeem <- redeem_df_ks %>% filter(statistic > 0.35) %>% pull(variant)
  biased_variants_mgatk <- mgatk_df_ks %>% filter(statistic > 0.35) %>% pull(variant)
  
  # Function to plot positions of alt alleles on molecule
  make_plot <- function(var_df, var_levels, color, title) {
    n_vars = length(unique(var_df$variant))
    pP <- ggplot(var_df %>%  mutate(variant2 = factor(as.character(variant), levels = var_levels)), 
                 aes( x = rel_position)) +
      geom_histogram(color = color, fill = color, bins = 30) + 
      facet_wrap(~variant2, scales = "free_y") +
      ggtitle(paste0(sample_id, "- ",title," top variants")) + theme_bw() +
      labs(x = "Position on molecule", y = "abundance of alt allele") 
    
    ggsave(pP, file = paste0("../plots_split/", sample_id, "_",title,".png"),
           width = (n_vars)^(1/3)*4, height = (n_vars)^(1/3)*2, dpi = 200) # semi arbitrary way of making plots ok
    title
  }
  
  pRlow <- make_plot(
    sum_stats_redeem %>% filter( !(variant %in% biased_variants_redeem)),
    top_hetero_variants_redeem, "dodgerblue3", "Redeem_KSretained")
  
  pRhigh <- make_plot(
    sum_stats_redeem %>% filter(variant %in% biased_variants_redeem),
    top_hetero_variants_redeem, "firebrick", "Redeem_KSremoved")
  
  pMlow <- make_plot(
    sum_stats_mgatk %>% filter( !(variant %in% biased_variants_mgatk)),
    top_hetero_variants_mgatk, "dodgerblue", "mgatk_KSretained")
  
  pMhigh <- make_plot(
    sum_stats_mgatk %>% filter( variant %in% biased_variants_mgatk),
    top_hetero_variants_mgatk, "red", "mgatk_KSremoved")
  
}

# Do one per donor
process_one_sample("Young1.T1.BMMC")
process_one_sample("Youn2.BMMC")
process_one_sample("Old1.BMMC")
process_one_sample("Old2.BMMC")

# See if donor-specific biases are accounted for in weighted jaccard

# Global: import mutation weights for jaccard statistic
data(CellPCT)
V.weight <- data.frame(weight=1-CellPCT$muRate)
V.weight$Variants <- paste0("Variants",gsub("_","",CellPCT$Variant))
V.weight %>% 
  filter(Variants %in% c("Variants2AG", "Variants3TA", "Variants4CT"))


