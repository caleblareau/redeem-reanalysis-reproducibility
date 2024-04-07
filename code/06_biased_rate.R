library(BuenColors)
source("00_functions.R")

# update this based on your computer
# https://figshare.com/articles/dataset/ReDeeM_raw_mutation_calling/24418966/1
# Directory should have these files from the published paper
base_dir <- "/Users/lareauc/Downloads/redeem-downloaded/mito_data_redeem/"  

# Main Function that processes summary counts
process_one_sample <- function(sample_id){
  
  # import redeem df
  redeemV_data <- fread(paste0(base_dir,sample_id,".Consensus.final/RawGenotypes.Sensitive.StrandBalance"))
  
  # import and select top 60 variants per variant caller
  sens_call_df <- fread(paste0("../data/redeem-sensitive-calls/",sample_id,".Consensus.final_S.tsv"))
  top_hetero_variants_redeem <- sens_call_df %>% filter(totalVAF < 0.8) %>%
    arrange(desc(CellN)) %>% pull(Variants)
  
  ## mgatk calls 
  mgatk_df_in <- fread(paste0("../data/mgatk-varstats/D_",sample_id,"_hg38_v20-mtMask.variant_stats.tsv.gz"))
  mgatk_df <- mgatk_df_in %>% filter(strand_correlation > 0.65 & n_cells_conf_detected >= 3 & (log10(vmr) > -2) & mean < 0.8)
  top_hetero_variants_mgatk <- mgatk_df %>% arrange(desc(n_cells_conf_detected)) %>% 
    mutate(redeem_formatted = paste0(position, "_", gsub(">", "_", nucleotide))) %>%
    pull(redeem_formatted) 
  
  # Now filter the redeemV data frame and annotate molecule position
  sum_stats_redeem <- redeemV_data %>% filter(V4 %in% top_hetero_variants_redeem) %>% make_position_df()
  sum_stats_mgatk <- redeemV_data %>% filter(V4 %in% top_hetero_variants_mgatk) %>% make_position_df()
  
  # do KS tests on all variants
  redeem_df_ks <- sum_stats_redeem %>% make_ks_test_df
  mgatk_df_ks <- sum_stats_mgatk %>% make_ks_test_df
  
  # Create a set of variants with bias per caller
  output_df <-  data.frame(
    variant = c(redeem_df_ks$variant, mgatk_df_ks$variant),
    statistic = c(redeem_df_ks$statistic, mgatk_df_ks$statistic),
    method = c(rep("redeem", length(redeem_df_ks$variant)), rep("mgatk", length(mgatk_df_ks$variant)))
  ) %>% mutate(biased = statistic > 0.35)  %>%
    group_by(method) %>% summarize(Percent_biased = mean(biased)*100)
  output_df %>% mutate(sample_id)
}

# Do one per donor
variant_filter_stats <- rbind(
  process_one_sample("Young1.T1.BMMC"),
  process_one_sample("Youn2.BMMC"),
  process_one_sample("Old1.BMMC"),
  process_one_sample("Old2.BMMC")
)

