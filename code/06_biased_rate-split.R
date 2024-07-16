library(BuenColors)
source("00_functions.R")

# Main Function that processes summary counts
process_one_sample <- function(sample_id){
  
  # import redeem df
  redeemV_data <- fread(paste0(base_dir,sample_id,".Consensus.final/RawGenotypes.Sensitive.StrandBalance"))
  
  # import and select variants per caller
  sens_call_df <- fread(paste0("../data/redeem-sensitive-calls/",sample_id,".Consensus.final_S.tsv"))
  top_hetero_variants_redeem <- sens_call_df %>% filter(totalVAF < 0.8) %>%
    arrange(desc(CellN)) %>% pull(Variants)
  
  sens_df_mut <- sens_call_df %>% filter(totalVAF < 0.8) %>%
    mutate(avg_support_per_pos_cell = TotalVcount/CellN)  %>%
    mutate(LMHC = avg_support_per_pos_cell < 1.3 & CellN > 20)
  lmhc_vars <- sens_df_mut %>% filter(LMHC) %>% pull(Variants)
  other_redeem_vars <- sens_df_mut %>% filter(!LMHC) %>% pull(Variants)
  
  ## mgatk calls 
  mgatk_df_in <- fread(paste0("../data/mgatk-varstats/D_",sample_id,"_hg38_v20-mtMask.variant_stats.tsv.gz")) %>%
    mutate(redeem_formatted = paste0(position, "_", gsub(">", "_", nucleotide)))
  mgatk_df <- mgatk_df_in %>% filter(strand_correlation > 0.65 & n_cells_conf_detected >= 3 & (log10(vmr) > -2) & mean < 0.9)
  very_high <- mgatk_df %>% filter(mean> 0.4) %>% pull(redeem_formatted)
  
  top_hetero_variants_mgatk <- mgatk_df %>% arrange(desc(n_cells_conf_detected)) %>% 
    pull(redeem_formatted) 
  
  homoplasmic_vars <- mgatk_df_in %>% arrange(desc(n_cells_conf_detected)) %>%
    filter(mean > 0.9) %>%
    mutate(redeem_formatted = paste0(position, "_", gsub(">", "_", nucleotide))) %>%
    pull(redeem_formatted) 
  
  # Now filter the redeemV data frame and annotate molecule position
  sum_stats_redeem_lmhc <- redeemV_data %>% filter(V4 %in% lmhc_vars) %>% 
    filter(!c(V4 %in% bad_vars_cw) & !c(V4 %in% very_high)) %>%
    make_position_df()
  sum_stats_redeem_others <- redeemV_data %>% filter(V4 %in% other_redeem_vars) %>% 
    filter(!c(V4 %in% bad_vars_cw) & !c(V4 %in% very_high)) %>%
    make_position_df()
  sum_stats_mgatk <- redeemV_data %>% filter(V4 %in% top_hetero_variants_mgatk) %>%
    filter(!c(V4 %in% bad_vars_cw) & !c(V4 %in% very_high)) %>%
    make_position_df()
  sum_stats_homoplasmic <- redeemV_data %>% filter(V4 %in% homoplasmic_vars) %>%
    filter(!c(V4 %in% bad_vars_cw) & !c(V4 %in% very_high)) %>%
    make_position_df()
  
  
  # Look at the proportion of each variant type
  var_lmhc <- unique(sum_stats_redeem_lmhc$variant)
  df_lmhc <- data.frame(
    table(substr(var_lmhc, nchar(var_lmhc)-2, nchar(var_lmhc)))
  ) %>% mutate(lmhc_pct = round(Freq/length(var_lmhc)*100, 1))
  
  var_mgatk <- unique(sum_stats_mgatk$variant)
  df_mgatk <- data.frame(
    table(substr(var_mgatk, nchar(var_mgatk)-2, nchar(var_mgatk)))
  ) %>% mutate(mgatk_pct = round(Freq/length(var_mgatk)*100, 1))
  
  # Merge the two data frames
  mdf <- merge(df_lmhc[,c(1,3)], df_mgatk[,c(1,3)], by = "Var1", all = TRUE)
  mdf[is.na(mdf)] <-  0
  mdf$transition <- mdf$Var1 %in% c("A_G", "C_T", "T_C", "G_A")
  mdf$variant <- gsub("_", ">", mdf$Var1)
  ggplot(mdf %>% reshape2::melt(id.vars = c("Var1", "variant", "transition")), 
         aes(x = variable, y = value, fill = transition)) +
    geom_bar(stat = "identity")
  print(paste0("# mgatk vars: ", length(var_mgatk)))
  print(paste0("# lmhc vars: ", length(var_lmhc)))
  
  # do KS tests on all variants
  redeem_others_df_ks <- sum_stats_redeem_others %>% make_ks_test_df
  redeem_LMHC_df_ks <- sum_stats_redeem_lmhc %>% make_ks_test_df
  mgatk_df_ks <- sum_stats_mgatk %>% make_ks_test_df
  homoplasmic_df_ks <- sum_stats_homoplasmic %>% make_ks_test_df
  
  # Create a set of variants with bias per caller
  output_df <-  data.frame(
    variant = c(redeem_others_df_ks$variant,redeem_LMHC_df_ks$variant, mgatk_df_ks$variant, homoplasmic_df_ks$variant),
    statistic = c(redeem_others_df_ks$statistic,redeem_LMHC_df_ks$statistic, mgatk_df_ks$statistic, homoplasmic_df_ks$statistic),
    method = c(rep("redeem_others", length(redeem_others_df_ks$variant)),rep("redeem_LMHC", length(redeem_LMHC_df_ks$variant)),
               rep("mgatk", length(mgatk_df_ks$variant)),  rep("homoplasmic", length(homoplasmic_df_ks$variant)))
  ) %>% mutate(biased = statistic > 0.35)
  output_df %>% mutate(sample_id)
}

# Do one per donor
variant_filter_stats <- rbind(
  process_one_sample("Young1.T1.BMMC"),
  process_one_sample("Youn2.BMMC"),
  process_one_sample("Old1.BMMC"),
  process_one_sample("Old2.BMMC")
)

# fix typo 
variant_filter_stats$sample_id <- gsub("Youn2.BMMC", "Young2.BMMC", variant_filter_stats$sample_id)

po <- ggplot(variant_filter_stats, aes(x = statistic, fill = method)) + 
  geom_histogram(bins = 21) + facet_grid(method ~ sample_id, scales = "free_y") +
  geom_vline(xintercept = 0.35, linetype = 2) +
  pretty_plot(fontsize = 8) + geom_hline(yintercept = 0) +
  labs(x = "KS statistic", y = "Abundance of variants") +
  scale_fill_manual(values = c("black", "dodgerblue3", "firebrick", "grey")) +
  theme(legend.position = "none")
cowplot::ggsave2(po, file = "../final_plots/ks_test_statistics.pdf", width = 7, height = 4)

variant_filter_stats %>%
  group_by(method, sample_id) %>%
  summarize(count = n(), n_biased = sum(biased), perc_b = round(mean(biased)*100, 1))

variant_filter_stats %>%
  filter(method == "mgatk" & biased) %>% arrange(desc(variant))

write.table(variant_filter_stats, file = "../output/variant_statisticss_KS_all.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)