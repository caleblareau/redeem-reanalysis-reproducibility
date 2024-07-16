library(BuenColors)
source("../../code/00_functions.R")

# Convert ATAC and RNA barcodes
atac_bc <- fread("../../data/multiome_barcodes/WhiteList_10X_Multiome.ATAC.gz", header = FALSE)[[1]]
rna_bc <- fread("../../data/multiome_barcodes/WhiteList_10X_Multiome.RNA.gz", header = FALSE)[[1]]
names(rna_bc) <- atac_bc

# import redeem df
what = "aged"
what = "young"

if(what == "young"){
  redeemV_data <- fread("../../../redeem-downloaded/mito_data_redeem/GSE261078_downloaded/GSM8133910_Extended_Young_BMMC.mtDNA.RawGenotypes.Total.StrandBalance.txt")
  donor_id <- fread("../../../redeem-downloaded/mito_data_redeem/GSE261078_downloaded/GSM8133912_Extended_Young_BMMC.Hash_calls.tsv")
} else {
  redeemV_data <- fread("../../../redeem-downloaded/mito_data_redeem/GSE261078_downloaded/GSM8133916_Extended_Aged_BMMC_HSPC.mtDNA.RawGenotypes.Total.StrandBalance.txt.gz")
  donor_id <- fread("../../../redeem-downloaded/mito_data_redeem/GSE261078_downloaded/GSM8133918_Extended_Aged_BMMC_HSPC.Hash_calls.tsv")
}

call_vec <- donor_id$call; names(call_vec) <- donor_id[["cell"]]
redeemV_data$donor <- call_vec[rna_bc[as.character(redeemV_data$V2)]]

# find all positions
all_frag_pos <- data.frame(
  first = str_split_fixed(redeemV_data[["V1"]], "_", 3)[,c(2)] %>% as.numeric(),
  last = str_split_fixed(redeemV_data[["V1"]], "_", 3)[,c(3)] %>% as.numeric(),
  pos = redeemV_data$V3, variant = redeemV_data$V4, donor = redeemV_data[["donor"]],
  cb = redeemV_data[["V2"]]
) %>% mutate(length = abs(last - first))  %>%
  filter(!(variant %in% c(bad_vars_cw, "3107_N_C", "3107_N_T", "3107_N_A", "3107_N_G")))

# Pull data out
get_edgeFC <- function(ggploto){
  df <-  ggplot_build(ggploto)$data[[1]]
  mean(df[c(1:2, 50,51),"y"])/mean(df[c(4:48),"y"])
}
get_edge_perc <- function(ggploto){
  df <-  ggplot_build(ggploto)$data[[1]]
  sum(df[c(1:2, 50,51),"y"])/sum(df[,"y"])*100
}

# Split by donor
donor_list <- split(all_frag_pos, all_frag_pos$donor)
#  one_donor_df <- donor_list[[1]]
lapply(1:length(donor_list), function(i){
  one_donor_df <- donor_list[[i]]
  j <- i+2
  
  summary_df <- one_donor_df %>%
    group_by(variant, cb) %>% summarize(obs_count = n()) %>%
    ungroup() %>% group_by(variant) %>% summarize(TotalVcount = sum(obs_count), CellN = sum(obs_count >= 0)) %>%
    arrange(desc(CellN)) %>% mutate(pct_cells_positive = CellN/max(CellN)) %>%
    mutate(avg_support_per_pos_cell = TotalVcount/CellN)
  lmic <- summary_df %>% 
    filter(avg_support_per_pos_cell <= 1.3 & CellN >= 20) %>% pull(variant)
  p1 <- one_donor_df %>% filter((variant %in% c( lmic))) %>%
    mutate(length = abs(last - first)) %>% 
    mutate(dist_to_edge = pmin(abs(last - pos), abs(first - pos)),
           pos_on_molecule = ifelse(abs(last - pos) < abs(first - pos), length - dist_to_edge, dist_to_edge)) %>%
    mutate(rel_position = pos_on_molecule/length) %>%
    ggplot(aes(x = rel_position)) + 
    geom_histogram(bins = 51, fill = "firebrick") + scale_y_continuous(expand = c(0,0)) + 
    labs(x = "position on molecule", y = "count of reference mismatch") + pretty_plot(fontsize = 6) + L_border()
  p1
  cowplot::ggsave2(p1, file = paste0("../plots/",what,"_",as.character(j), ".pdf"), width = 1.5, height = 1.5)
  
  data.frame(
    what = c( "redeem"),
    value = round(c(get_edgeFC(p1)), 2),
    percent_edge = round(c(get_edge_perc(p1)), 2)
  )
  
})

