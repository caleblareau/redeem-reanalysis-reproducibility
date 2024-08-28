library(data.table)
library(dplyr)
library(redeemR)
require(dplyr)
require(Matrix.utils)
library(ggtree)
library(ape)
library(phangorn)
library(ggtreeExtra)
source("00_functions.R")

# example
id = "Young1.T1.BMMC"
compute_trees_v1_v2 <- function(id){
  print(id)
  WD <-  paste0(base_dir, id, ".Consensus.final/") # CL specific working directory
  
  # Import the old way
  redeemR<-Create_redeemR(redeemR.read(path=WD,thr="S",Processed=F,rdsname = "/new.VariantsGTSummary.RDS"))
  redeemR@HomoVariants <- (redeemR@V.fitered %>% filter(totalVAF > 0.35) %>% pull(Variants)) 
  
  # New method -- trimming
  VariantsGTSummary_trim5 <- redeemR.read.trim(path=WD, "S", Processed=F,
                                               rdsname="/VariantsGTSummary.S.trim5_binom.RDS",edge_trim=5)
  
  redeemR_trim5_binom<-Create_redeemR_model(VariantsGTSummary_trim5)
  redeemR_trim5_binom<- clean_redeem(redeemR_trim5_binom,fdr = 0.05)
  redeemR_trim5_binom<-clean_redeem_removehot(redeemR_trim5_binom)
  
  # Append the variants that Chen manually filters out already as homoplasmic
  # so they aren't in the matrix
  redeemR@HomoVariants <- c(redeemR@HomoVariants, bad_vars_cw2)
  
  # Create matrix using redeem functions
  redeemR <- Make_matrix(redeemR)
  
  ## Filter low coverage cells and intersect with trim5
  BadCells<-subset(redeemR@CellMeta,meanCov<10)$Cell
  keep_cells <- rownames(redeemR@Cts.Mtx.bi)[!(rownames(redeemR@Cts.Mtx.bi) %in% BadCells)]
  keep_cells2 <- intersect(keep_cells,
                           intersect(rownames(redeemR_trim5_binom@Cts.Mtx.bi), rownames(redeemR@Cts.Mtx.bi)))
  length(keep_cells2)
  Cts.Mtx.bin1 <- redeemR@Cts.Mtx[keep_cells2,]
  Cts.Mtx.bin1[Cts.Mtx.bin1>=1]<-1
  
  Cts.Mtx.bin2 <- redeemR_trim5_binom@Cts.Mtx[keep_cells2,!(Variants_to_underscores( colnames(redeemR_trim5_binom@Cts.Mtx.bi)) %in% redeemR@HomoVariants)]
  Cts.Mtx.bin2[Cts.Mtx.bin2>=1]<-1
  
  mean(rowSums(Cts.Mtx.bin1)==1)
  mean(rowSums(Cts.Mtx.bin2)==1)
  
  # function to make the weights for their jaccard thingy
  make_weight_df <- function(Cts.bin){
    # Pull out weights fromt the data object; update with N/As with 1 as described in manuscript
    weightdf<-data.frame(Variants=colnames(Cts.bin)) %>% merge(.,V.weight,by="Variants",all.x = T,sort = F) 
    weight <- weightdf %>% pull(weight)
    weight[is.na(weight)] <- 1
    weight
  }
  
  # do the weights
  weight1 <- make_weight_df(Cts.Mtx.bin1)
  weight2 <- make_weight_df(Cts.Mtx.bin2)
  
  # Compute weighted jaccard distance metric in the paper
  w_j1 <- quick_w_jaccard_cl(Cts.Mtx.bin1,weight1) %>% data.matrix()
  w_j2 <- quick_w_jaccard_cl(Cts.Mtx.bin2,weight2) %>% data.matrix()

  # Make phylogenetic tree
  phylo1 <- nj(w_j1)
  phylo2 <- nj(w_j2)
  
  save(phylo1, phylo2,  file = paste0("../output/trees_version/trees_v1_v2_", id, ".rda"))
  
}

compute_trees_v1_v2("Youn2.BMMC")
compute_trees_v1_v2("Young1.T1.BMMC")
compute_trees_v1_v2("Old1.BMMC")
compute_trees_v1_v2("Old2.BMMC")



make2trees_redeem_version <- function(id, w_t_2){
  
  load(paste0("../output/trees_version/trees_v1_v2_", id, ".rda"))
  
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
  

  
  cowplot::ggsave2(t1f, file = paste0("../final_plots/tree_version/tree_redeemV1_", id, ".pdf"),
                   width = 3, height = 3)
  cowplot::ggsave2(t2f, file = paste0("../final_plots/tree_version/tree_redeemV2_", id, ".pdf"),
                   width = 3, height = 3)
}

make2trees_redeem_version("Youn2.BMMC", 10)
make2trees_redeem_version("Young1.T1.BMMC", 5)

make2trees_redeem_version("Old1.BMMC", 10)
make2trees_redeem_version("Old2.BMMC", 7)

