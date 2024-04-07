library(BuenColors)
source("00_functions.R")

# Example
sample_id = "Youn2.BMMC"
sample_id = "Young1.T1.BMMC"
sample_id = "Old2.BMMC"

# Global vars
base_dir <- "/Users/lareauc/Downloads/redeem-downloaded/mito_data_redeem/"
bad_vars_cw <- c("310_T_C","3109_T_C","309_C_T", "5764_C_T") # variants that we agree are bad should always be filtered

# Import data
redeemV_data <- fread(paste0(base_dir,sample_id,".Consensus.final/RawGenotypes.Sensitive.StrandBalance"))
mgatk_df_in <- fread(paste0("../data/mgatk-varstats/D_",sample_id,"_hg38_v20-mtMask.variant_stats.tsv.gz"))
mgatk_df <- mgatk_df_in %>% filter(strand_correlation > 0.65 & n_cells_conf_detected >= 3 & (log10(vmr) > -2) & mean < 0.8)
mgatk_vars <- mgatk_df %>% arrange(desc(n_cells_conf_detected)) %>% 
  mutate(redeem_formatted = paste0(position, "_", gsub(">", "_", nucleotide))) %>%
  pull(redeem_formatted) 

# Load scatters
redeem_df <- fread(paste0("../data/redeem-sensitive-calls/",sample_id,".Consensus.final_S.tsv"))

# Define thresholds of suspicous variants
threshold_avg <- 1.3
threshold_cells <- 30

sens_df_annotated <- redeem_df %>%
  filter(!c(Variants %in% bad_vars_cw)) %>%
  mutate(avg_support_per_pos_cell = TotalVcount/CellN) %>%
  mutate(color = case_when(Variants %in% mgatk_vars ~ "mgatkcall",
                           CellN > 30 & avg_support_per_pos_cell < threshold_avg  ~ "redeemProblem",
                           TRUE ~ "redeemOK"))

ggplot(sens_df_annotated %>% arrange(desc(color)), aes(x = CellN, y = avg_support_per_pos_cell, label = Variants, color = color)) +
  geom_point() + scale_x_log10() + scale_y_log10() +
  geom_hline(yintercept = threshold, linetype = 2)

sens_df_annotated %>% filter(CellN > 30 & color == "mgatkcall" & avg_support_per_pos_cell < 1.3)
sens_df_annotated %>% filter(CellN > 30 & avg_support_per_pos_cell < 1.3) %>% dim()

redeemV_data %>%
  filter(V4 %in% (sens_df_annotated %>% filter(color == "redeemProblem") %>% pull(Variants))) %>%
  make_position_df() %>%
  ggplot(aes(x = rel_position)) + 
  geom_histogram(bins = 101) +
  ggtitle('problematic') + labs(x = "length of molecule") + theme_bw()

#Make a full profile
redeemV_data %>%
  make_position_df() %>%
  ggplot(aes(x = rel_position)) + 
  geom_histogram(bins = 101) 
