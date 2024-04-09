library(Seurat)
source("00_functions.R")
list.files(paste0(base_dir, "/../seurat_data_redeem"))

so <- readRDS(paste0(base_dir, "/../seurat_data_redeem/Young1.All.T1.Seurat.RDS"))
d1_meta <- so@meta.data[,c("STD.CellType","STD_Cat","ClonalGroup","Sample")]

so2 <- readRDS(paste0(base_dir, "/../seurat_data_redeem/Young2.All.Seurat.RDS"))
d2_meta <- so2@meta.data[,c("STD.CellType","STD_Cat","ClonalGroup","Sample")]

soO1 <- readRDS(paste0(base_dir, "/../seurat_data_redeem/Old1.BMMC_HSPC.Seurat.RDS"))
soO1_meta <- soO1@meta.data[,c("STD.CellType","STD_Cat","ClonalGroup","Sample")]

soO2 <- readRDS(paste0(base_dir, "/../seurat_data_redeem/Old2.BMMC_HSPC.Seurat.RDS"))
soO2_meta <- soO2@meta.data[,c("STD.CellType","STD_Cat","ClonalGroup","Sample")]

write.table(d1_meta, file = "../output/Young1_seurat_meta.tsv",
            row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)

write.table(d2_meta, file = "../output/Young2_seurat_meta.tsv",
            row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)

write.table(soO1_meta, file = "../output/Old1_seurat_meta.tsv",
            row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)

write.table(soO2_meta, file = "../output/Old2_seurat_meta.tsv",
            row.names = TRUE, col.names = FALSE, sep = "\t", quote = FALSE)


