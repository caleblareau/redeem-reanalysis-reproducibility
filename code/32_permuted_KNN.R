# Most of this code is from  https://github.com/chenweng1991/redeem_robustness_reproducibility/blob/master/notebooks/notebook_share/06.KNN_celltype_origin.ipynb
# Modifications are required as the processed data used for the real analysis 
# were not provided (i.e. whatever "~/data/Young1_BMMC_HSPC_Redeem_filter1.rds" is)
# 

library(redeemR)
library(ggplot2)
library(dplyr)
library(reshape2)
library(Matrix)
library(amap)
library(gridExtra)
source("00_functions.R")


# From https://github.com/chenweng1991/redeemR/issues/5
# How to combine multiple redeem objects once imported
combined_redeemR<-function(ob.list,thr="S"){
  ## Combine GTsummary(Add suffix at the end of cell names)
  mutated_list <- lapply(seq_along(ob.list), function(i) {
    x <- ob.list[[i]]@GTsummary.filtered
    suffix <- paste0("_", i)
    x <- mutate(x, Cell = paste0(Cell, suffix))
    return(x)
  })
  combined.GTsummary <- do.call(rbind, mutated_list)
  ## Combine depth
  pos_cov<-ob.list[[1]]@DepthSummary[[1]][,1,drop=F]
  pos_cov$meanCov<-lapply(ob.list,function(x){x@DepthSummary[[1]][,2]}) %>% do.call(cbind,.) %>% rowSums
  ## Combine depth at cell level(Add suffix at the end of cell names)
  mutated_list <- lapply(seq_along(ob.list), function(i) {
    x <- ob.list[[i]]@DepthSummary[[2]]
    suffix <- paste0("_", i)
    x <- mutate(x, V1 = paste0(V1, suffix),Sample=names(ob.list)[i])
    return(x)
  })
  cell_cov<-do.call(rbind, mutated_list)
  combined.DepthSummary<-list(pos_cov,cell_cov)
  names(combined.DepthSummary)<-c("Pos.MeanCov","Cell.MeanCov")
  ## Add attributes
  attr(combined.GTsummary,"depth")<-combined.DepthSummary
  attr(combined.GTsummary, "thr")<-thr
  combined_redeemR<-Create_redeemR(combined.GTsummary)
  return(combined_redeemR)
}

# Import meta data and translate the two
rna2atac <- fread("../data/multiome_barcodes/WhiteList_10X_Multiome.ATAC.gz", header = FALSE)[[1]]
names(rna2atac) <- fread("../data/multiome_barcodes/WhiteList_10X_Multiome.RNA.gz", header = FALSE)[[1]]
meta_df <- fread("../output/Young1_seurat_meta.tsv") 
ct <- as.character(meta_df$V2); names(ct) <- paste0(unname(rna2atac[substr( meta_df$V1,1, 16)]), "_", substr(meta_df$V1,18, 18))

# Filter 2
# Import Young 1 BMMC and HPC
id = "Young1.T1.BMMC"
bmmc_VariantsGTSummary_trim5 <- redeemR.read.trim(path=paste0(base_dir, id, ".Consensus.final/"), "S", Processed=F,
                                                  rdsname="/VariantsGTSummary.S.trim5_binom.RDS",edge_trim=5)

bmmc_redeemR_trim5_binom<-Create_redeemR_model(bmmc_VariantsGTSummary_trim5)
bmmc_redeemR_trim5_binom<- clean_redeem(bmmc_redeemR_trim5_binom,fdr = 0.05)
bmmc_redeemR_trim5_binom<-clean_redeem_removehot(bmmc_redeemR_trim5_binom)

id = "Young1.T1.HPC"
hpc_VariantsGTSummary_trim5 <- redeemR.read.trim(path=paste0(base_dir, id, ".Consensus.final/"), "S", Processed=F,
                                                 rdsname="/VariantsGTSummary.S.trim5_binom.RDS",edge_trim=5)

hpc_redeemR_trim5_binom<-Create_redeemR_model(hpc_VariantsGTSummary_trim5)
hpc_redeemR_trim5_binom<- clean_redeem(hpc_redeemR_trim5_binom,fdr = 0.05)
hpc_redeemR_trim5_binom<-clean_redeem_removehot(hpc_redeemR_trim5_binom)

Young1_BMMC_HSPC_redeem.sensitive <- combined_redeemR(list(bmmc_redeemR_trim5_binom,hpc_redeemR_trim5_binom), "S")

# Import data based on reproducibility notebook
Young1_BMMC_HSPC_redeem.sensitive@CellMeta$STD.CellType <- unname(ct[Young1_BMMC_HSPC_redeem.sensitive@CellMeta$Cell])
Young1_BMMC_HSPC_redeem.sensitive <- Make_matrix(Young1_BMMC_HSPC_redeem.sensitive, )
keep_cells <- intersect(
  Young1_BMMC_HSPC_redeem.sensitive@CellMeta[complete.cases(Young1_BMMC_HSPC_redeem.sensitive@CellMeta), "Cell"][[1]],
  rownames(Young1_BMMC_HSPC_redeem.sensitive@Cts.Mtx.bi)
)


## Knn analysis function from Chen
draw_heatmap<-function(redeem, filterV, cell_meta_df, shuffle=1000,cut=0.05,k=14){
  M<-redeem@Cts.Mtx.bi[keep_cells,!colnames(redeem@Cts.Mtx.bi) %in% paste("Variants",gsub("_","",filterV),sep="")]
  M<-M[rowSums(M)>0,]
  M<-M[,colSums(M)>0]    
  print(paste("Number of variants to compute distance:",ncol(M)))
  D<-as.matrix(BinaryDist(M,method = "Jaccard"))
  
  ## Generate celltype.mtx which is a cell vs celltype matrix based on knn
  celltype.mtx<-c()
  Rnames<-row.names(D)
  CellDic<-cell_meta_df[,"STD.CellType"]; CellDic[CellDic == "Refined.HSC"] <- "HSC"
  names(CellDic)<-cell_meta_df$Cell  
  Cnames<-unique(sort(CellDic))
  
  #sapply to determine how many of each celltype is around the current cell
  celltype.mtx <- sapply(1:ncol(D), function(i){
    current.cell<-colnames(D)[i]
    name.ordered<-D[,i] %>% .[order(.)] %>% names
    x<-setdiff(name.ordered,current.cell)[1:k]
    count_vec <- as.numeric(table(CellDic[x])[Cnames])
    count_vec[is.na(count_vec)] <- 0
    count_vec
  })
  celltype.mtx <- t(celltype.mtx)
  row.names(celltype.mtx)<-Rnames
  colnames(celltype.mtx)<-Cnames
  celltype.df<-as.data.frame(celltype.mtx)
  celltype.df$Seed=CellDic[Rnames]
  
  #Compute aggregated cell type number for each neighbourhood
  NeibourComp_mean<-celltype.df %>% group_by(Seed) %>% dplyr::summarise(MKP=mean(MKP),CMP=mean(CMP),MPP=mean(MPP),HSC=mean(HSC),LMPP=mean(LMPP),Mono=mean(Mono),GMP=mean(GMP),MDP=mean(MDP), cDC=mean(cDC),MEP=mean(MEP),EryP=mean(EryP),ProB=mean(ProB),CLP=mean(CLP),pDC=mean(pDC),B=mean(B),NK=mean(NK),CD8=mean(CD8),CD4=mean(CD4)) 
  cell_order<-c("Mono","MDP","GMP","cDC","EryP","MEP","MKP","HSC","MPP","CMP","pDC","LMPP","CLP","ProB","B","NK","CD8","CD4")
  NeibourComp_mean$Seed<-factor(NeibourComp_mean$Seed,levels=cell_order)
  
  #Normalization
  m.normed<-NeibourComp_mean[,2:ncol(NeibourComp_mean)] %>% apply(.,2,function(x){x/median(x)}) %>% cbind(NeibourComp_mean[,1],.)
  m.normed$Seed<-paste("S.",m.normed[,1],sep="")
  m.normed<-m.normed[!m.normed$Seed %in% c("S.HSC","S.MPP","S.MKP","S.LMPP","S.CMP"),] %>% dplyr::select(-HSC,-MPP,-LMPP,-MKP,-CMP)
  m.normed.obs<-m.normed
  hmod<-m.normed.obs %>% tibble::remove_rownames() %>% tibble::column_to_rownames("Seed") %>% Dist(.,method = "correlation") %>% hclust 
  
  ## Do reshuffle
  celltype.df.reshuffle<-celltype.df
  m.normed.reshuffle.list<-list()
  for(i in 1:shuffle){
    celltype.df.reshuffle$Seed<-sample(celltype.df.reshuffle$Seed)
    NeibourComp_mean<-celltype.df.reshuffle %>% group_by(Seed) %>% dplyr::summarise(MKP=mean(MKP),CMP=mean(CMP),MPP=mean(MPP),HSC=mean(HSC),LMPP=mean(LMPP),Mono=mean(Mono),GMP=mean(GMP),MDP=mean(MDP), cDC=mean(cDC),MEP=mean(MEP),EryP=mean(EryP),ProB=mean(ProB),CLP=mean(CLP),pDC=mean(pDC),B=mean(B),NK=mean(NK),CD8=mean(CD8),CD4=mean(CD4)) 
    NeibourComp_mean$Seed<-factor(NeibourComp_mean$Seed,levels=cell_order)
    m.normed<-NeibourComp_mean[,2:ncol(NeibourComp_mean)] %>% apply(.,2,function(x){x/median(x)}) %>% cbind(NeibourComp_mean[,1],.)
    m.normed$Seed<-paste("S.",m.normed[,1],sep="")
    #Compute m.normed
    m.normed<-m.normed[!m.normed$Seed %in% c("S.HSC","S.MPP","S.MKP","S.LMPP","S.CMP"),] %>% dplyr::select(-HSC,-MPP,-LMPP,-MKP,-CMP)
    m.normed.reshuffle.list<-c(m.normed.reshuffle.list,list(m.normed))
  }
  
  #Compute p value
  mtx<-m.normed.obs[,-1]
  pvalue.mtx<-matrix(NA, nrow = nrow(mtx), ncol = ncol(mtx))
  for (i in 1:nrow(mtx)){
    for(j in 1:ncol(mtx)){
      result_vector <- sapply(m.normed.reshuffle.list, function(x) x[i, (j+1)])
      p<-(length(result_vector)-(sum(mtx[i,j]>result_vector)))/length(result_vector)
      pvalue.mtx[i,j]<-p
    }
  }
  colnames(pvalue.mtx)<-colnames(mtx)
  pvalue.mtx<-as.data.frame(pvalue.mtx)
  pvalue.mtx$Seed<-m.normed.obs$Seed
  pvalue.mtx.m<-melt(pvalue.mtx)
  pvalue.mtx.m$ID<-paste(pvalue.mtx.m$Seed,pvalue.mtx.m$variable,sep="_")     
  pvalue.mtx.m$value[pvalue.mtx.m$value>cut]<-NA
  pvalue.mtx.m$value<-as.character(pvalue.mtx.m$value)
  pvalue.mtx.m$value[pvalue.mtx.m$value=="0"]<-"***"
  #
  datatoplot<-melt(m.normed.obs)
  ## prepare datatoplot
  seed.l<-c('S.GMP','S.pDC','S.MEP','S.EryP','S.MDP','S.cDC','S.Mono','S.CD8','S.CD4','S.B','S.NK','S.CLP','S.ProB')
  variable.l<-c('ProB','CLP','NK','B','CD4','CD8','Mono','cDC','MDP','EryP','MEP','pDC','GMP')                                
  
  datatoplot$Seed<-factor(datatoplot$Seed,levels=seed.l)
  datatoplot$variable<-factor(datatoplot$variable,levels=variable.l)
  datatoplot$ID<-paste(datatoplot$Seed,datatoplot$variable,sep="_")
  datatoplot<-merge(datatoplot,pvalue.mtx.m[,c("ID","value")],by="ID")
  p.main<-ggplot(datatoplot)+aes(variable,Seed,fill=log(value.x),label=value.y)+geom_tile()+scale_fill_gradient2(high="darkred",low="steelblue",mid="white",midpoint=0,na.value = "steelblue",limits=c(-1,1),oob=scales::squish)+
    theme_classic()+theme(axis.text = element_text(size=16,color="black"))
  
  return(p.main)
}

# New function to compute lineage-non-informative variants
# key workflow from C. Weng but put in function to facilitate permutation
determine_non_informative_variants <- function(){
  
  # Define a data frame of meta data since missing data in Redeem object
  cell_meta_df <- data.frame(Young1_BMMC_HSPC_redeem.sensitive@CellMeta)
  rownames(cell_meta_df) <- cell_meta_df$Cell
  cell_meta_df <- cell_meta_df[keep_cells,]
  
  
  set.seed(5)
  cell_meta_df$STD.CellType<- sample(cell_meta_df$STD.CellType)
  
  
  lineage_vector <- dplyr::recode(cell_meta_df$STD.CellType,
                                  MK="MKP", 
                                  EryP ="ME",MEP="ME",
                                  HSC="HSCMPP", MPP="HSCMPP",LMPP="HSCMPP", CMP="HSCMPP", Refined.HSC  = "HSCMPP",
                                  MDP="Mye", GMP="Mye",cDC="Mye", Mono="Mye",
                                  ProB="Lym", CLP="Lym", pDC="Lym",B="Lym", NK="Lym", CD8="Lym", CD4="Lym")
  cell_meta_df$Lineage <- lineage_vector
  
  ## Prepare variants summary
  LV<-Young1_BMMC_HSPC_redeem.sensitive@GTsummary.filtered %>% 
    merge(.,cell_meta_df[,c("Cell","Lineage")]) 
  LV <- LV[complete.cases(LV),]
  
  ## Make LV.summarise
  LV.summarise<-LV %>% group_by(Variants) %>% dplyr::count(Lineage) %>% dcast(Variants~Lineage)
  row.names(LV.summarise)<-LV.summarise$Variants
  LV.summarise<-LV.summarise[,-1]
  LV.summarise[is.na(LV.summarise)]<-0
  LV.summarise$sum<-rowSums(LV.summarise)
  
  ## A function to test lineage distribution
  VariantLinDistri.compute<-function(x,alt="greater"){
    ps<-c()
    Ess<-c()
    for(i in 1:5){
      p<-binom.test(x[i],x[6],Exp[i],alternative = alt)$p.value
      Es<-(x[i]/x[6])/Exp[i]
      ps<-c(ps,p)
      Ess<-c(Ess,Es)
    }
    return(list(ps,Ess)) 
  }
  ## Calculate expected lineage proportion
  Exp<-colSums(LV.summarise[,-6])/sum(colSums(LV.summarise[,-6]))
  
  # Compute pvalues
  LV.summarise.raw<-apply(LV.summarise,1,VariantLinDistri.compute)
  LV.summarise.enrich.pvalues<- lapply(LV.summarise.raw,function(x){x[[1]]}) %>% do.call(rbind,.) %>% as.data.frame
  LV.summarise.enrich.effectSizes<-lapply(LV.summarise.raw,function(x){x[[2]]}) %>% do.call(rbind,.)  %>% as.data.frame
  
  non_info_variants<-LV.summarise.enrich.effectSizes[which((LV.summarise.enrich.effectSizes$Mye<=1.5 | LV.summarise.enrich.pvalues$Mye>=0.05) &
                                                             (LV.summarise.enrich.effectSizes$Lym<=1.5 | LV.summarise.enrich.pvalues$Lym>=0.05) &
                                                             (LV.summarise.enrich.effectSizes$ME<=1.5  | LV.summarise.enrich.pvalues$ME>=0.05)  &
                                                             (LV.summarise.enrich.effectSizes$MKP<=1.5  | LV.summarise.enrich.pvalues$MKP>=0.05)),] %>% row.names 
  
  sum_stats <- data.frame(n_non_info = length(non_info_variants), n_tested = length(rownames(LV.summarise.enrich.effectSizes))) %>%
    mutate(n_informative = n_tested - n_non_info)
  print(sum_stats)
  return(list(cell_meta_df, non_info_variants))
}

# Make function calls that do the heavy lifting
permuted_niv <- determine_non_informative_variants()
permuted_plot <- draw_heatmap(Young1_BMMC_HSPC_redeem.sensitive, filterV = permuted_niv[[2]], cell_meta_df = permuted_niv[[1]])

cowplot::ggsave2(permuted_plot+ggtitle("PERMUTED DATA with lineage-informative filtering") + geom_text(), 
                 file = "../final_plots/permuted_knn.pdf", width = 8, height = 7)
