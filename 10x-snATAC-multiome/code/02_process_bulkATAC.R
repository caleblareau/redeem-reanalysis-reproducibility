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
cdf <- fread("../data/SRR7245894_bulkATAC.MD.tsv.gz") %>% head(200000)
cdf$toParse <- gsub("MD:Z:", "", cdf$V3)

# remove indels and Ns for simplicity
cdf <- cdf[!grepl("\\^", cdf$toParse) & !grepl("I", cdf$toParse) & !grepl("N", cdf$toParse),]

# split forward and reverse and call base quality
cdf$mate <- bamFlagAsBitMatrix(cdf$V1)[,"isSecondMateRead"] 
cdf$strand <- bamFlagAsBitMatrix(cdf$V1)[,"isMinusStrand"] 

cdf_r1 <- cdf %>% filter(cdf$strand == 1 & mate == 0) %>% process_one_df
cdf_r2 <- cdf %>% filter(cdf$strand == 1 & mate == 1) %>% process_one_df
cdf_f1 <- cdf %>% filter(cdf$strand == 0 & mate == 0) %>% process_one_df
cdf_f2 <- cdf %>% filter(cdf$strand == 0 & mate == 1) %>% process_one_df

save(cdf_r1, cdf_r2, cdf_f1, cdf_f2, file = "../output/summary_position_bias_bulkATAC.rda")


# Now visualize data
load("../output/summary_position_bias_bulkATAC.rda")

p_r2 <- ggplot(cdf_r2 %>% filter(bq >= 30), aes(x = idx)) + 
  geom_histogram(bins = 38) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(1,25,50)) +
  labs(x = "position on read", y = "# alternate alleles")

p_f2 <- ggplot(cdf_f2 %>% filter(bq >= 30), aes(x = idx)) + 
  geom_histogram(bins = 38) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = c(1,25,50)) +
  labs(x = "position on read", y = "# alternate alleles")

p_r1 <- ggplot(cdf_r1 %>% filter(bq >= 30), aes(x = idx)) + 
  geom_histogram(bins = 38) +  
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(breaks = c(1,25,50)) +
  labs(x = "position on read", y = "# alternate alleles")
p_f1 <- ggplot(cdf_f1 %>% filter(bq >= 30), aes(x = idx)) + 
  geom_histogram(bins = 38) +  
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(expand = c(0,0)) + 
  scale_x_continuous(breaks = c(1,25,50)) +
  labs(x = "position on read", y = "# alternate alleles")

cowplot::ggsave2(cowplot::plot_grid(p_r1, p_f1,
                                    p_r2, p_f2), file = "../output/bulkATAC_full.pdf", 
                 width = 3.5, height = 2.5)
