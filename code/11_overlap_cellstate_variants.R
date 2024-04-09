library(data.table)
library(dplyr)
library(redeemR)
library(matrixStats)

source("00_functions.R")
atac_bc <- fread("../data/multiome_barcodes/WhiteList_10X_Multiome.ATAC.gz", header = FALSE)[[1]]
rna_bc <- fread("../data/multiome_barcodes/WhiteList_10X_Multiome.RNA.gz", header = FALSE)[[1]]
names(rna_bc) <- atac_bc


# example
id = "Young1.T1.BMMC"
meta_df <- fread("../output/Young1_seurat_meta.tsv",
                 col.names = c("barcode", "STD.CellType","STD_Cat","ClonalGroup","Sample")) %>% data.frame()

analyze_cell_state_clade_impact <- function(id, meta_df){
  WD <-  paste0(base_dir, id, ".Consensus.final/") # CL specific working directory
  
  # Import data based on reproducibility notebook
  redeemR<-Create_redeemR(redeemR.read(path=WD,thr="S",Processed=F,rdsname = "/new.VariantsGTSummary.RDS"))
  redeemR@HomoVariants <- (redeemR@V.fitered %>% filter(totalVAF > 0.35) %>% pull(Variants)) # Variants in Old 2 that are at ~40%
  
  # Append the variants that Chen manually filters out already as homoplasmic
  # so they aren't in the matrix
  redeemR@HomoVariants <- c(redeemR@HomoVariants, bad_vars_cw)
  
  # Identify potential problematic variants from the KS test
  redeemV_data <- fread(paste0(WD, "/RawGenotypes.Sensitive.StrandBalance"))
  redeem_called_vars <- redeemR@V.fitered$Variants
  sum_stats_redeem <- redeemV_data %>% filter(V4 %in% redeem_called_vars) %>% 
    make_position_df() %>% make_ks_test_df()
  suspected_bad <- sum_stats_redeem %>% filter(statistic > 0.35) %>% pull(variant)
  suspected_bad_formatted <- paste0("Variants", gsub("_", "", suspected_bad))
  
  # Create matrix using redeem functions
  redeemR <- Make_matrix(redeemR)
  bad_vars_cw %in% colSums(redeemR@Cts.Mtx.bi) 
  
  ## Filter low coverage cells
  BadCells<-subset(redeemR@CellMeta,meanCov<10)$Cell
  keep_cells <- !(rownames(redeemR@Cts.Mtx.bi) %in% BadCells)
  Cts.Mtx.bin <- redeemR@Cts.Mtx.bi[keep_cells,]
  
  # now filter out the position biased variants
  boo_keep <- !(colnames(Cts.Mtx.bin) %in% suspected_bad_formatted)
  sum(boo_keep)
  
  # Assess which idx in seurat meta file matches here
  idxs <- sort(unique(substr(meta_df$barcode, 18, 18)))
  idx_go <- sapply(idxs, function(idx_1){
    sum(paste0(rna_bc[rownames(Cts.Mtx.bin)], "-", idx_1) %in% meta_df$barcode)
  }) %>% which.max()
  
  updated_names <- paste0(rna_bc[rownames(Cts.Mtx.bin)], "-", idx_go)
  barcodes_use <- intersect(updated_names, meta_df$barcode)
  stopifnot(length(barcodes_use) > 1000)
  
  # Create a map back
  atac_bc_x <- paste0(atac_bc, "-", idx_go)
  names(atac_bc_x) <- paste0(rna_bc, "-", idx_go)
  
  # Now subset to cells jointly analyzable
  rownames(Cts.Mtx.bin) <- paste0(rownames(Cts.Mtx.bin), "-", idx_go)
  Cts.Mtx.bin.ss <- Cts.Mtx.bin[unname(atac_bc_x[barcodes_use]),]
  cell_state_df <- meta_df; rownames(cell_state_df) <- cell_state_df$barcode
  cell_state_df <- cell_state_df[barcodes_use,]
  clusters_ordered <- cell_state_df$STD.CellType
  Cat_ordered <- cell_state_df$STD_Cat
  ClonalGroup_ordered <- as.character(cell_state_df$ClonalGroup)
  
  # Do the cell type biased first
  var_meta_df <- sum_stats_redeem
  var_meta_df$variantF <-  paste0("Variants", gsub("_", "", var_meta_df$variant))
  common_vars <- intersect(colnames(Cts.Mtx.bin.ss)[colSums(Cts.Mtx.bin.ss) > 1], var_meta_df$variantF)
  var_meta_df <- var_meta_df %>% filter(variantF %in% common_vars)
  Cts.Mtx.bin.ss <- Cts.Mtx.bin.ss[,common_vars]
  
  # Now loop and compute a pvalue
  chisq_pvalue_df <- lapply(1:dim(var_meta_df)[1], function(i){ # 
    var_df <- data.frame(
      cluster = clusters_ordered,
      cat = Cat_ordered,
      clone = ClonalGroup_ordered,
      variant = Cts.Mtx.bin.ss[,var_meta_df$variantF[i]] , stringsAsFactors = TRUE
    )
    
    # fix instances where
    # test is undefined b/c no variants (rare)
    res.cs1 <- chisq.test(var_df$variant, var_df$cluster)
    res.cs2 <- chisq.test(var_df$variant, var_df$cat)
    res.cs3 <- chisq.test(var_df$variant,var_df$clone)
    cluster_pvalue = res.cs1[[3]]
    lineage_pvalue = (res.cs2)[[3]]
    clone_pvalue = (res.cs3)[[3]]
    
    # Also pull out the top clade
    top_clade <- sort(table(ClonalGroup_ordered[ Cts.Mtx.bin.ss[,var_meta_df$variantF[i]] > 0 ]), decreasing = TRUE) %>% head(1)
    
    data.frame(
      variant = var_meta_df$variantF[i],
      cluster_pvalue,
      lineage_pvalue,
      clone_pvalue,
      n_vars = sum(var_df$variant),
      top_clade = names(top_clade),
      prop_vars_in_clade = as.numeric(unname(top_clade)/sum(var_df$variant)),
      prop_clade_with_var = as.numeric(unname(top_clade)/(sum(names(top_clade) == ClonalGroup_ordered, na.rm=TRUE)))
    )
  }) %>% rbindlist()
  
  mdf <- merge(chisq_pvalue_df, var_meta_df, by.x = "variant", by.y = "variantF")
  mdf %>% arrange((clone_pvalue)) 
}

col.namesV = c("barcode", "STD.CellType","STD_Cat","ClonalGroup","Sample")

y1df <- analyze_cell_state_clade_impact("Young1.T1.BMMC", fread("../output/Young1_seurat_meta.tsv", col.names = col.namesV ) %>% data.frame())
y2df <- analyze_cell_state_clade_impact("Youn2.BMMC", fread("../output/Young2_seurat_meta.tsv", col.names = col.namesV ) %>% data.frame())
o1df <- analyze_cell_state_clade_impact("Old1.BMMC", fread("../output/Old1_seurat_meta.tsv", col.names = col.namesV ) %>% data.frame())
o2df <- analyze_cell_state_clade_impact("Old2.BMMC", fread("../output/Old2_seurat_meta.tsv", col.names = col.namesV ) %>% data.frame())

write.table(y1df, file = "../output/Young1_variant_statistics.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(y2df, file = "../output/Young2_variant_statistics.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(o1df, file = "../output/Old1_variant_statistics.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(o2df, file = "../output/Old2_variant_statistics.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

qplot(y2df$statistic, y = -log10(y2df$cluster_pvalue))

y2df %>% filter(statistic > 0.35 & clone_pvalue < 10^-10)

y2df %>% filter(statistic > 0.35 & cluster_pvalue < 10^-3)

ggplot(y2df, aes(x = statistic, y = prop_clade_with_var)) +
  geom_point()
