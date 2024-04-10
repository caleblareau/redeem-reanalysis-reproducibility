library(data.table)
library(dplyr)
library(stringr)
library(BuenColors)
library(Rsamtools)
 

parse_MD_row <- function(row_one) {
  all_matches <- str_extract_all(row_one$toParse, "[0-9]+[ACGT]")[[1]] # tailing numbers we don't care about
  mismatches <- as.integer(sub("[ACGT]", "", all_matches))
  starting_indices <- as.integer(gregexpr("[0-9]+", all_matches)[[1]]) 
  idx <- starting_indices + cumsum(mismatches)
  bq <- sapply(idx, function(i){
    as.numeric(charToRaw(substr(row_one$V2, i, i))) - 33
  })
  return(data.frame(idx, bq))
}

process_one_df <- function(df){
  lapply(1:dim(df)[1], function(i){
    parse_MD_row(df[i,])
  }) %>% rbindlist() %>% data.frame()
}

# Pull in data
cdf <- fread("chr22_bwa.tsv")
cdf$toParse <- gsub("MD:Z:", "", cdf$V3)

# remove indels and Ns for simplicity
cdf <- cdf[!grepl("\\^", cdf$toParse) & !grepl("I", cdf$toParse) & !grepl("N", cdf$toParse),]

# split forward and reverse and call base quality
cdf$strand <- bamFlagAsBitMatrix(cdf$V1)[,"isMinusStrand"] 
cdf_f <- cdf %>% filter(cdf$strand == 1) %>% process_one_df
cdf_r <- cdf %>% filter(cdf$strand == 0) %>% process_one_df

p1 <- ggplot(cdf_f %>% filter(bq >= 30), aes(x = idx)) + 
  geom_histogram(bins = 50) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(1,25,50)) +
  labs(x = "position on read", y = "# alternate alleles")

p2 <- ggplot(cdf_r %>% filter(bq >= 30), aes(x = idx)) + 
  geom_histogram(bins = 50) +  
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(breaks = c(1,25,50)) +
  labs(x = "position on read", y = "# alternate alleles")

cowplot::ggsave2(p1, file = "scATAC_positiveStrand.pdf", 
                 width = 2, height = 1.5)
cowplot::ggsave2(p2, file = "scATAC_negativeStrand.pdf", 
                 width = 2, height = 1.5)
