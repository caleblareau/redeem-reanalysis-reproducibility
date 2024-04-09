library(data.table)
library(dplyr)
library(qvalue)
library(ggbeeswarm)
col.namesV = c("barcode", "STD.CellType","STD_Cat","ClonalGroup","Sample")

process_donor <- function(id){
  # Import clonal clades
  sm <-fread(paste0("../output/",id,"_seurat_meta.tsv"), col.names = col.namesV)
  uniq_clonal_groups <- sort(unique(sm$ClonalGroup)); ugc <- ugc[!is.na(ugc)]
  
  # Import variant associations
  stat_df <- fread(paste0("../output/",id,"_variant_statistics.tsv"))
  stat_df$qvalue_clade <- qvalue(stat_df$clone_pvalue)[["qvalues"]]
  stat_df$qvalue_cluster <- qvalue(stat_df$cluster_pvalue)[["qvalues"]]
  
  rdf1 <- stat_df %>% filter(qvalue_cluster < 0.05) %>%
    summarize(n = n(), n_biased = sum(statistic > 0.35)) %>%
    mutate(what = "celltype_cluster", who = id)
  
  rdf2 <- stat_df %>% filter(qvalue_clade < 0.05) %>%
    summarize(n = n(), n_biased = sum(statistic > 0.35)) %>%
    mutate(what = "clade", who = id)
  rbind(rdf1, rdf2)
}
lapply(c("Young1", "Young2", "Old1", "Old2"), process_donor) %>%
  rbindlist() %>% mutate(perc = n_biased/n*100) %>%
  ggplot( aes(x = what, y = perc, color = who)) +
  geom_quasirandom(width = 0.1) +
  labs(x = "", y = "% of informative variants with position bias") + 
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(limits = c(0,100)) +
  geom_hline(yintercept = c(0, 100), linetype = 2) + 
  scale_color_manual(values = jdb_palette("corona")[c(4,3,1,2)]) +
  theme(legend.position = "none") -> p1
cowplot::ggsave2(p1, file = "../final_plots/informative_variant_filtered.pdf", width = 1.5, height = 1.7)


                         