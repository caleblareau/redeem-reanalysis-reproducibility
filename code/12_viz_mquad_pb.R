library(BuenColors)
source("00_functions.R")

# Example
sample_id = "Old2.BMMC"

process_one_mquad_calls <- function(sample_id){
  
  # Import data
  redeemV_data <- fread(paste0(base_dir,sample_id,".Consensus.final/RawGenotypes.Sensitive.StrandBalance"))
  mgatk_df_in <- fread(paste0("../data/mgatk-varstats/D_",sample_id,"_hg38_v20-mtMask.variant_stats.tsv.gz"))
  mgatk_df <- mgatk_df_in %>% filter(strand_correlation > 0.65 & n_cells_conf_detected >= 3 & (log10(vmr) > -2) & mean < 0.9) %>%
    mutate(redeem_formatted = paste0(position, "_", gsub(">", "_", nucleotide)))
  mgatk_vars <- mgatk_df %>% arrange(desc(n_cells_conf_detected)) %>% 
    pull(redeem_formatted) 
  very_high <- mgatk_df %>% filter(mean> 0.4) %>% pull(redeem_formatted)
  
  # Definie homoplasmic
  homoplasmic_vars <- mgatk_df_in %>% arrange(desc(n_cells_conf_detected)) %>%
    filter(mean > 0.9) %>%
    mutate(redeem_formatted = paste0(position, "_", gsub(">", "_", nucleotide))) %>%
    pull(redeem_formatted) 
  
  # Import mquad data
  bic <- fread(paste0("../../redeem-downloaded/large_for_mquad/",sample_id,"/BIC_params.csv"))
  variant_list <- fread(paste0("../../redeem-downloaded/large_for_mquad/",sample_id,".variants.out.tsv"), header = FALSE)[[1]]
  snp_idx <- bic %>% filter(PASS_KP) %>%
    pull(variant_name) %>% gsub("SNP", "", .) %>% as.character() %>% as.numeric()
  mquad_vars <- variant_list[snp_idx] %>% Variants_to_underscores
  mquad_vars <- mquad_vars[!c(mquad_vars %in% c(homoplasmic_vars, very_high, bad_vars_cw))]
  print(length(mquad_vars))
  
  # Load scatters
  redeem_df <- fread(paste0("../data/redeem-sensitive-calls/",sample_id,".Consensus.final_S.tsv"))
  
  # Define thresholds of suspicous variants
  threshold_avg <- 1.3
  threshold_cells <- 20
  
  sens_df_annotated <- redeem_df %>%
    filter(!c(Variants %in% bad_vars_cw) & !c(Variants %in% very_high)) %>%
    mutate(avg_support_per_pos_cell = TotalVcount/CellN) %>%
    mutate(color = case_when(Variants %in% mquad_vars ~ "amquadcall",
                             CellN > threshold_cells & avg_support_per_pos_cell < threshold_avg  ~ "LMHC",
                             TRUE ~ "OK"))
  
  pA1 <- ggplot(sens_df_annotated %>% arrange((color)), aes(x = CellN, y = avg_support_per_pos_cell, label = Variants, color = color)) +
    geom_point(size = 0.5) + scale_x_log10() + scale_y_log10() +
    geom_hline(yintercept = threshold_avg, linetype = 2) +
    pretty_plot(fontsize = 6) + labs(x = "# cells with mutation", y = "mean # molecules / + cell") +
    scale_color_manual(values = c("purple3", "firebrick","lightgrey")) +
    theme(legend.position = "none") + L_border()
  
  # Make Mquad pileup plot
  pMQ <- redeemV_data %>%
    filter(V4 %in% mquad_vars) %>%
    make_position_df() %>%
    ggplot(aes(x = rel_position)) + 
    geom_histogram(bins = 51, fill = "purple3") + scale_y_continuous(expand = c(0,0)) + 
    labs(x = "position on molecule", y = "count of reference mismatch") + pretty_plot(fontsize = 6) + L_border()
  
  pMQ
  cowplot::ggsave2(cowplot::plot_grid(pA1, pMQ, nrow = 2),
                   width = 1.5, height = 3, file = paste0("../final_plots/", sample_id, "_mQUAD.pdf"))
  
  # Pull data out
  get_edgeFC <- function(ggploto){
    df <-  ggplot_build(ggploto)$data[[1]]
    mean(df[c(1:3, 49,50,51),"y"])/mean(df[c(4:48),"y"])
  }
  get_edge_perc <- function(ggploto){
    df <-  ggplot_build(ggploto)$data[[1]]
    sum(df[c(1:3, 49,50,51),"y"])/sum(df[,"y"])*100
  }
  data.frame(
    what = c("mQuad" ),
    value = round(c(get_edgeFC(pMQ)), 2),
    percent_edge = round(c(get_edge_perc(pMQ)), 2)
  ) 
  
}

process_one_mquad_calls( "Young1.T1.BMMC") # 
process_one_mquad_calls( "Old1.BMMC") # 
process_one_mquad_calls( "Youn2.BMMC") # 
process_one_mquad_calls( "Old2.BMMC") # 

