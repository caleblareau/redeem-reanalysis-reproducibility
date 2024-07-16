library(BuenColors)
source("00_functions.R")

# Example
sample_id = "Youn2.BMMC"

process_one_bar_prop <- function(sample_id){
  
  # import and select variants per caller
  sens_call_df <- fread(paste0("../data/redeem-sensitive-calls/",sample_id,".Consensus.final_S.tsv"))
  top_hetero_variants_redeem <- sens_call_df %>% filter(totalVAF < 0.8) %>%
    arrange(desc(CellN)) %>% pull(Variants)
  
  sens_df_mut <- sens_call_df %>% filter(totalVAF < 0.8) %>%
    mutate(avg_support_per_pos_cell = TotalVcount/CellN)  %>%
    mutate(LMHC = avg_support_per_pos_cell < 1.3 & CellN > 20)
  
  ## mgatk calls 
  mgatk_df_in <- fread(paste0("../data/mgatk-varstats/D_",sample_id,"_hg38_v20-mtMask.variant_stats.tsv.gz")) %>%
    mutate(redeem_formatted = paste0(position, "_", gsub(">", "_", nucleotide)))
  mgatk_df <- mgatk_df_in %>% filter(strand_correlation > 0.65 & n_cells_conf_detected >= 3 & (log10(vmr) > -2) & mean < 0.9)
  very_high <- mgatk_df %>% filter(mean> 0.4) %>% pull(redeem_formatted)
  
  # Redeem data
  redeemV_data <- fread(paste0(base_dir,sample_id,".Consensus.final/RawGenotypes.Sensitive.StrandBalance")) %>%
    filter(V4 %in% sens_df_mut$Variants) %>%
    mutate(barcode_variant = paste0(V2, "_", V4)) 
  pos_all_redeem <- make_position_df(redeemV_data)
  redeemV_data$rel_position <- pos_all_redeem$rel_position
  
  redeemV_bv_count_change <- function(countidx){
    
    # filter by positions on molecule
    what_changes <- redeemV_data %>% 
      filter(rel_position >= 0.95 | rel_position <= 0.05) %>%
      mutate(mut = paste0(V5, "_", V6)) %>%
      group_by(V2, V4, mut) %>% summarize(count = n()) %>%
      filter(count == countidx) %>%
      pull(mut)
    
    # What frequency of nucleotide changes occur
    data.frame(table(what_changes)) %>%
      mutate(prop = Freq/sum(Freq)) %>% mutate(countidx)
  }
  
  redeemV_bv_count_change_min <- function(countidx){
    
    # filter by positions on molecule
    what_changes <- redeemV_data %>% 
      filter(rel_position >= 0.95 | rel_position <= 0.05) %>%
      mutate(mut = paste0(V5, "_", V6)) %>%
      group_by(V2, V4, mut) %>% summarize(count = n()) %>%
      filter(count >= countidx) %>%
      pull(mut)
    
    # What frequency of nucleotide changes occur
    data.frame(table(what_changes)) %>%
      mutate(prop = Freq/sum(Freq)) %>% mutate(countidx)
  }
  
  
  lapply(1:2, redeemV_bv_count_change) %>%
    rbindlist() %>% data.frame() -> odf
  
  odf <- rbind(odf, redeemV_bv_count_change_min(3))
  trav <- odf %>% filter(!(what_changes %in% c("C_T", "A_G", "G_A", "T_C")))
  print(trav %>% group_by(countidx) %>% summarize(sum(prop)))
  p1 <- ggplot(trav, aes(x = countidx, y = prop*100, fill = what_changes)) + 
    geom_bar(stat = "identity", color = "black", width = 0.8)  + scale_fill_manual(values = jdb_palette("corona")) + 
    scale_y_continuous(expand = c(0,0), breaks = c(0,10,20,30), limits = c(0, 35)) + pretty_plot(fontsize = 7) + L_border() +
    labs(x = "# eUMIs per variant + cell", y = "% of all substitutions", color = "") +
    scale_x_continuous(breaks = c(1:6)) + theme(legend.position = "none")
  p1
  cowplot::ggsave2(p1, file = paste0("../final_plots/mut_sig_", sample_id, "_123.pdf"), 
                   width = 1.4, height = 1.8)
  
}

process_one_bar_prop( "Young1.T1.BMMC") 
process_one_bar_prop( "Old1.BMMC") 
process_one_bar_prop( "Youn2.BMMC")  
process_one_bar_prop( "Old2.BMMC")  

