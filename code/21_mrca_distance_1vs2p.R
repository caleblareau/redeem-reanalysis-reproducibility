library(castor)
library(dplyr)
library(data.table)
library(BuenColors)

if(FALSE){
  process_tree_MRCA <- function(t1, t2){
    n_nodes <- t1$Nnode
    
    lapply(1:100, function(j){
      lapply(1:1000, function(i){
        set.seed(i+j*50000)
        random3 <- sort(sample(1:n_nodes, size = 3))
        eg <- expand.grid(random3, random3) %>%
          filter(Var1 < Var2)
        d1 <- get_pairwise_distances(t1, A = eg[,1], B= eg[,2])
        d2 <- get_pairwise_distances(t2, A = eg[,1], B= eg[,2])
        data.frame(idx = i,
                   perm = j,
                   match = which.min(d1) == which.min(d2))
      }) %>% rbindlist() %>% data.frame() -> odf
      odf
    }) %>% rbindlist() %>% data.frame() ->  full_df
    full_df
  }
  
  process_donor <- function(id){
    load(paste0("../output/trees_partition/trees_all_1_2p_",id,".rda"))
    mrca_all_vs_2p <- process_tree_MRCA(phylo1, phylo2)
    mrca_all_vs_1 <- process_tree_MRCA(phylo1, phylo3)
    rbind(mrca_all_vs_2p %>% mutate(comparison = "all_vs_2p"),
          mrca_all_vs_1 %>% mutate(comparison = "all_vs_1")    )
  }
  young1 <- process_donor("Young1.T1.BMMC")
  young2 <- process_donor("Youn2.BMMC")
  old1 <- process_donor("Old1.BMMC")
  old2 <- process_donor("Old2.BMMC")
  
  save(old1, old2, young1, young2, file = "../output/trees/tree_MRCA_stats_split1vs2p.rda")
  
}

##### 
# load and plot
####

load("../output/trees/tree_MRCA_stats_split1vs2p.rda")

rbind(
  old1 %>% mutate(donor = "zAged1"),
  old2 %>% mutate(donor = "Aged2"),
  young1 %>% mutate(donor = "zzYoung1"),
  young2 %>% mutate(donor = "zYoung2")
) %>%
  group_by(donor, comparison, idx) %>%
  summarize(prop = mean(match)) -> stats_df
mean_df <- stats_df %>% group_by(donor, comparison ) %>% summarize(prop = mean(prop))
mean_df
stats_df %>%
  ggplot(aes(x = donor, y= prop*100, color = comparison)) + 
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  coord_cartesian(ylim = c(25, 100)) +
  scale_color_manual(values = c("grey3","grey")) +
  geom_hline(yintercept = c(33), color = "firebrick", linetype = 2 )+ 
  geom_hline(yintercept = c(100), color = "dodgerblue3", linetype = 2 )+ 
  coord_flip() +
  pretty_plot(fontsize = 8) + labs(y = "% MRCA trio concordance") + L_border() +
  theme(legend.position = "none") -> px

px
cowplot::ggsave2(px ,file = "../final_plots/mrca_rates_split.pdf", width = 2.5, height = 2.5)
