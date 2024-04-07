library(BuenColors)
source("00_functions.R")

process_one_sample_percent_mutated <- function(id_one){
  n_positions_mutated <- fread(paste0("../data/redeem-sensitive-calls/",id_one,".Consensus.final_S.tsv")) %>%
    mutate(position = as.numeric(str_extract(.$Variants, pattern = "[0-9]+"))) %>%
    pull(position) %>% unique() %>% length() 
  data.frame(id_one, n_positions_mutated)
}


rbind(
  process_one_sample_percent_mutated("Young1.T1.BMMC"),
  process_one_sample_percent_mutated("Youn2.BMMC"),
  process_one_sample_percent_mutated("Old1.BMMC"),
  process_one_sample_percent_mutated("Old2.BMMC")
) %>% mutate(percent = n_positions_mutated/16569*100)
