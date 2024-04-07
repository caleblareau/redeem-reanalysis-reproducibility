library(data.table)
library(dplyr)
library(redeemR)
library(stringr)
library(tidyr)
library(purrr)
library(BuenColors)


# update this based on your computer
# https://figshare.com/articles/dataset/ReDeeM_raw_mutation_calling/24418966/1
# Directory should have files from the published paper downloaded from figshare
base_dir <- "/Users/lareauc/Dropbox/redeem/redeem-downloaded/mito_data_redeem/"  

# variants that we agree are bad / already filtered
bad_vars_cw <- c("310_T_C","3109_T_C","309_C_T", "5764_C_T") 


## Helper functions 
## loading packages

# Function that takes the Redeem V output dataframe
# and appends where on the molecule the alt allele is observed
# relative to the edges of the alignment
make_position_df <- function(in_df){
  data.frame(
    first = str_split_fixed(in_df[["V1"]], "_", 3)[,c(2)] %>% as.numeric(),
    last = str_split_fixed(in_df[["V1"]], "_", 3)[,c(3)] %>% as.numeric(),
    pos = in_df$V3, variant = in_df$V4
  ) %>% mutate(length = abs(last - first)) %>% 
    mutate(dist_to_edge = pmin(abs(last - pos), abs(first - pos)),
           pos_on_molecule = ifelse(abs(last - pos) < abs(first - pos), length - dist_to_edge, dist_to_edge)) %>%
    mutate(rel_position = pos_on_molecule/length) 
}

# Function to compute KS test statistics 
# For all variants given a redeem molecule
# data frame 
make_ks_test_df <- function(var_df){
  var_df %>%
  nest(data = c(-variant)) %>%
    mutate(D_stat=map(data, ~(ks.test(.x$rel_position, "punif",0,1))),
           tidied = purrr::map(D_stat, broom::tidy)) %>%
    mutate(n = map_dbl(data, nrow)) %>%
    tidyr::unnest(tidied) %>%
    dplyr::select(variant, statistic, p.value, n)
}


