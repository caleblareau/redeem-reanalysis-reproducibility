library(BuenColors)
source("../code/00_functions.R")

# import redeem df
redeemV_data <- fread("GSM8113080_Crispr_Mouse_Batch1.mtDNA.RawGenotypes.Total.StrandBalance.txt.gz")

# find all fragments with issues
all_frag_pos <- data.frame(
  first = str_split_fixed(redeemV_data[["V1"]], "_", 3)[,c(2)] %>% as.numeric(),
  last = str_split_fixed(redeemV_data[["V1"]], "_", 3)[,c(3)] %>% as.numeric(),
  pos = redeemV_data$V3, variant = redeemV_data$V4
) %>% mutate(length = abs(last - first), idx = 1:n()) 

all_frag_pos %>% filter(((pos > first) & (pos > last)) | ((pos < first) & (pos < last))) %>%
  pull(idx) -> which_problem_idx

length(which_problem_idx) / dim(all_frag_pos)[1]



all_frag_pos %>% filter(((pos > first) & (pos > last)) | ((pos < first) & (pos < last))) %>%
  arrange(desc(first)) %>% filter(last < 9000 & first < 9000 ) %>% pull(idx)   -> selected_really_bad

redeemV_data[selected_really_bad,] %>% filter(V7 > 3) %>% head()

# Some variants are called outside of where the fragment can be determined, making this 
# data problematic in some aspects; making a note but moving on
redeemV_data %>% filter(V4 == "4287_A_T") %>% filter(grepl("5661|24231", V1))



# Nevertheless, try to identify mutations that are pervasive in different cells
n_cells <- length(unique(redeemV_data$V2))

# Define thresholds of suspicous variants
threshold_avg <- 1.3
threshold_cells <- 20

sum_stats_df <- redeemV_data %>% 
  group_by(V2, V4) %>% summarize(obs_count = n()) %>%
  ungroup() %>% group_by(V4) %>% summarize(TotalVcount = sum(obs_count), CellN = sum(obs_count >= 0)) %>%
  arrange(desc(CellN)) %>% mutate(pct_cells_positive = CellN/n_cells) %>%
  mutate(avg_support_per_pos_cell = TotalVcount/CellN) %>%
  mutate(color = case_when( CellN > threshold_cells & avg_support_per_pos_cell < threshold_avg  ~ "redeemProblem",
                            TRUE ~ "redeemOK"))

# Make similar plot
pa <- ggplot(sum_stats_df %>% arrange(desc(color)), aes(x = CellN, y = avg_support_per_pos_cell, color = color, label = V4)) +
  geom_point(size = 0.5) + scale_x_log10() + scale_y_log10() +
  geom_hline(yintercept = threshold_avg, linetype = 2) +
  pretty_plot(fontsize = 6) + labs(x = "# cells with mutation", y = "mean # molecules / + cell") +
  scale_color_manual(values = c( "lightgrey","firebrick")) +
  theme(legend.position = "none") + L_border()

pz <- redeemV_data %>%
  filter(V4 %in% (sum_stats_df %>% filter(color == "redeemProblem") %>% pull(V4))) %>%
  make_position_df() %>%
  filter(rel_position >= 0 & rel_position <= 1) %>%
  ggplot(aes(x = rel_position)) + 
  geom_histogram(bins = 51, fill = "firebrick") + scale_y_continuous(expand = c(0,0)) + 
  labs(x = "position on molecule", y = "count of alternate allele") + pretty_plot(fontsize = 6) + L_border()
pz
get_edgeFC_mouse <- function(ggploto){
  df <-  ggplot_build(ggploto)$data[[1]]
  mean(df[c(1,2, 50,51),"y"])/mean(df[c(3:49),"y"])
}
data.frame(
  what = c("redeem"),
  value = round(c(get_edgeFC_mouse(pz)), 2)
) 

cowplot::ggsave2(pa, 
                 width = 1.5, height = 1.5, file = paste0("../final_plots/mouse_redeem_plotA.pdf"))


cowplot::ggsave2(pz, 
                 width = 1.5, height = 1.5, file = paste0("../final_plots/mouse_redeem_plotZ.pdf"))

