library(data.table)
library(dplyr)
library(redeemR)
library(stringr)
library(BuenColors)
library(ggtree)
library(ape)
library(phangorn)
library(ggtreeExtra)
source('../code/00_functions.R')


make3trees <- function(id, w_t_2){
  
  load(paste0("../output/trees_partition/trees_all_1_2p_", id, ".rda"))
  
  t1 <- ggtree(phylo1,layout="fan", branch.length='none')
  ncell <- length(get_taxa_name(t1))
  clmdf <- data.frame(
    Cell = get_taxa_name(t1),
    idx = 1:(ncell),
    val = 1:(ncell)/(ncell)
  )
  
  # Make first tree
  t1f <- t1 + geom_fruit( 
    data=clmdf, 
    geom=geom_tile, 
    mapping=aes(y=Cell,x=2,fill=val), 
    pwidth=0.001, 
    width=10, 
    offset=0.05
  ) + scale_fill_gradientn(colors = jdb_palette("brewer_spectra")) + theme(legend.position = "none")
  
  # Make new tree after filtering
  t2 <- ggtree(phylo2,layout="fan", branch.length='none')
  
  t2f <- t2 + geom_fruit( 
    data=clmdf, 
    geom=geom_tile, 
    mapping=aes(y=Cell,x=2,fill=val), 
    pwidth=0.001, 
    width=w_t_2, 
    offset=0.05
  ) + scale_fill_gradientn(colors = jdb_palette("brewer_spectra")) + theme(legend.position = "none")
  
  
  # Make new tree after filtering
  t3 <- ggtree(phylo3,layout="fan", branch.length='none')
  
  t3f <- t3 + geom_fruit( 
    data=clmdf, 
    geom=geom_tile, 
    mapping=aes(y=Cell,x=2,fill=val), 
    pwidth=0.001, 
    width=15, 
    offset=0.05
  ) + scale_fill_gradientn(colors = jdb_palette("brewer_spectra")) + theme(legend.position = "none")
  
  
  cowplot::ggsave2(t1f, file = paste0("../3tree/tree_alll_", id, ".pdf"),
                   width = 3, height = 3)
  cowplot::ggsave2(t2f, file = paste0("../3tree/tree_2p_", id, ".pdf"),
                   width = 3, height = 3)
  cowplot::ggsave2(t3f, file = paste0("../3tree/tree_exactly1_", id, ".pdf"),
                   width = 3, height = 3)
  id
}

make3trees("Young1.T1.BMMC", 170)

make3trees("Youn2.BMMC", 170)
make3trees("Old2.BMMC", 150)
make3trees("Old1.BMMC", 170)
