library(BuenColors)
source("00_functions.R")

# Example
sample_id = "Youn2.BMMC"

process_me <- function(sample_id){
  
  redeemV_data <- fread(paste0(base_dir,sample_id,".Consensus.final/RawGenotypes.Sensitive.StrandBalance")) 
  redeemV_data$count <- redeemV_data %>% group_by(V2, V4) %>% mutate(count = n()) %>% pull(count)
  pos_df <- make_position_df(redeemV_data)
  redeemV_data$rel_position <- pos_df$rel_position
  redeemV_data$keep <-  (pos_df$rel_position >= 0.05 & pos_df$rel_position <= 0.95) | redeemV_data$count >= 3
  
  make_plot <- function(df, title){
  df %>%
    ggplot(aes(x = rel_position, fill = keep)) + 
    geom_histogram(bins = 51) + scale_y_continuous(expand = c(0,0)) + 
    labs(x = "position on molecule", y = "count of alternate allele") + pretty_plot(fontsize = 6) + L_border() +
    ggtitle(title) + 
    scale_fill_manual(values = c("TRUE" = "lightgrey", "FALSE" = "firebrick")) + 
    theme(legend.position = "none")
  }
  p1 <- make_plot(redeemV_data %>% filter(count == 1) , "==1")
  p2 <- make_plot(redeemV_data %>% filter(count == 2) , "==2")
  p3 <- make_plot(redeemV_data %>% filter(count > 2) , ">=3")
  
  cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, nrow = 1), 
                   width = 5, height = 1.5, file = paste0("../output/trees_partition/", sample_id, "_split_count_pileup_2color.pdf"))
  
}
process_me("Young1.T1.BMMC")
process_me("Youn2.BMMC")
