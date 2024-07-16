library(data.table)
library(dplyr)
library(Seurat)

# Downloaded seurat data from figshare:
# https://figshare.com/articles/dataset/Annotated_Seurat_objects/23290004/1

# data is coming from here: 
# https://github.com/chenweng1991/redeem_reproducibility/blob/master/Note-5%20HSC%20clonal%20structure%20and%20heterogeneity.ipynb
so1 <- readRDS("Young1.HSC.T1T2.Seurat.RDS")
head(so1@meta.data)
dim((so1))

table(substr(colnames(so1), 17,21))
saveRDS((so1@meta.data), file = "donor1_HSC.T1.T2.rds")

# https://github.com/chenweng1991/redeemR/issues/5
