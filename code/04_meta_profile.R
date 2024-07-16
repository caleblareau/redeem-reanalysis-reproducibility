library(BuenColors)
source("00_functions.R")

# Example
sample_id = "Old2.BMMC"

process_one <- function(sample_id){
  
  # Import data
  redeemV_data <- fread(paste0(base_dir,sample_id,".Consensus.final/RawGenotypes.Sensitive.StrandBalance"))
  mgatk_df_in <- fread(paste0("../data/mgatk-varstats/D_",sample_id,"_hg38_v20-mtMask.variant_stats.tsv.gz"))
  mgatk_df <- mgatk_df_in %>% filter(strand_correlation > 0.65 & n_cells_conf_detected >= 3 & (log10(vmr) > -2) & mean < 0.9) %>%
    mutate(redeem_formatted = paste0(position, "_", gsub(">", "_", nucleotide)))
  mgatk_vars <- mgatk_df %>% arrange(desc(n_cells_conf_detected)) %>% 
    pull(redeem_formatted) 
  very_high <- mgatk_df %>% filter(mean> 0.4) %>% pull(redeem_formatted)
  print("some variants were very high...")
  print(mgatk_df %>% filter(mean> 0.4))
  print(" ")
  
  homoplasmic_vars <- mgatk_df_in %>% arrange(desc(n_cells_conf_detected)) %>%
    filter(mean > 0.9) %>%
    mutate(redeem_formatted = paste0(position, "_", gsub(">", "_", nucleotide))) %>%
    pull(redeem_formatted) 
  
  # Load scatters
  redeem_df <- fread(paste0("../data/redeem-sensitive-calls/",sample_id,".Consensus.final_S.tsv"))
  
  # Define thresholds of suspicous variants
  threshold_avg <- 1.3
  threshold_cells <- 20
  
  sens_df_annotated <- redeem_df %>%
    filter(!c(Variants %in% bad_vars_cw) & !c(Variants %in% very_high)) %>%
    mutate(avg_support_per_pos_cell = TotalVcount/CellN) %>%
    mutate(color = case_when(Variants %in% mgatk_vars ~ "mgatkcall",
                             CellN > threshold_cells & avg_support_per_pos_cell < threshold_avg  ~ "redeemProblem",
                             TRUE ~ "redeemOK"))
  
  mgatk_overlap_problem <- sens_df_annotated %>% filter(CellN > threshold_cells & color == "mgatkcall" & avg_support_per_pos_cell < threshold_avg)
  all_problems <- sens_df_annotated %>% filter(CellN > threshold_cells & avg_support_per_pos_cell < threshold_avg) %>% dim()
  
  print("variants that mgatk called that we are claiming are problematic")
  print(mgatk_overlap_problem)
  
  print("other problematic variants")
  print(all_problems)
  
  
  p1 <- ggplot(sens_df_annotated %>% arrange(desc(color)), aes(x = CellN, y = avg_support_per_pos_cell, label = Variants, color = color)) +
    geom_point(size = 0.5) + scale_x_log10() + scale_y_log10() +
    geom_hline(yintercept = threshold_avg, linetype = 2) +
    pretty_plot(fontsize = 6) + labs(x = "# cells with mutation", y = "mean # molecules / + cell") +
    scale_color_manual(values = c("dodgerblue3", "lightgrey","firebrick")) +
    theme(legend.position = "none") + L_border()
  
  p2 <- redeemV_data %>%
    filter(V4 %in% (sens_df_annotated %>% filter(color == "redeemProblem") %>% pull(Variants))) %>%
    make_position_df() %>%
    ggplot(aes(x = rel_position)) + 
    geom_histogram(bins = 51, fill = "firebrick") + scale_y_continuous(expand = c(0,0)) + 
    labs(x = "position on molecule", y = "count of alternate allele") + pretty_plot(fontsize = 6) + L_border()
  
  p3 <- redeemV_data %>%
    filter(V4 %in% (sens_df_annotated %>% filter(color == "mgatkcall") %>% pull(Variants))) %>%
    make_position_df() %>%
    ggplot(aes(x = rel_position)) + 
    geom_histogram(bins = 51, fill = "dodgerblue3") + scale_y_continuous(expand = c(0,0)) + 
    labs(x = "position on molecule", y = "count of alternate allele") + pretty_plot(fontsize = 6) + L_border()
  
  #Make a profile of homoplasmic variants for control
  p4 <- redeemV_data %>%
    filter(V4 %in% homoplasmic_vars) %>%
    make_position_df() %>%
    ggplot(aes(x = rel_position)) + scale_y_continuous(expand = c(0,0)) + 
    geom_histogram(bins = 51, fill = "darkgrey") +
    labs(x = "position on molecule", y = "count of alternate allele") + pretty_plot(fontsize = 6)+ L_border()
  
  #cowplot::ggsave2(cowplot::plot_grid(p1, p4, p3, p2, nrow = 1), 
  #                 width = 5.5, height = 1.1, file = paste0("../final_plots/", sample_id, "_plots.pdf"))
  
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
   what = c("mgatk", "homoplasmic", "redeem_low"),
   value = round(c(get_edgeFC(p3), get_edgeFC(p4), get_edgeFC(p2)), 2),
   percent_edge = round(c(get_edge_perc(p3), get_edge_perc(p4), get_edge_perc(p2)), 2)
   
 ) 
  
}

process_one( "Young1.T1.BMMC") # 
process_one( "Old1.BMMC") # 
process_one( "Youn2.BMMC") # 
process_one( "Old2.BMMC") # 

