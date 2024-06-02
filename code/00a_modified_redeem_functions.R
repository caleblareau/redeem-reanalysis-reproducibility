source("00_functions.R")

# Additional parameter "remove edge" to remove the edge variants
# If false, creates a redeem object with ONLY edge
GTSummary_CL_edge_partition <-function(RawGenotypes,filterN=T, remove_edge = TRUE, ratio_based_filter = TRUE){
  # Make Depth dictionary
  data(ContextsDic)
  
  # pull out the position of the mutation on the allele
  RawGenotypes_cp <- RawGenotypes
  colnames(RawGenotypes_cp) <- paste0("V", 1:14)
  pos_df <- make_position_df(RawGenotypes_cp)
  
  # keep only things not at the edges and in 1 or two copies
  count <- RawGenotypes_cp %>% group_by(V2, V4) %>% mutate(count = n()) %>% pull(count)
  boo_away_from_edge_or_abundant <-  (pos_df$rel_position >= 0.05 & pos_df$rel_position <= 0.95) | count >= 3
  boo_away_from_edge_or_abundant_cw <-  (pos_df$dist_to_edge > 4) | count >= 3
  
  if(ratio_based_filter){
    boo_keep <- boo_away_from_edge_or_abundant
  } else {
    boo_keep <- boo_away_from_edge_or_abundant_cw
  }
  
  if(remove_edge){
    RawGenotypes <- RawGenotypes[boo_keep, ]
  } else {
    RawGenotypes <- RawGenotypes[!boo_keep, ]
  }
  
  
  Depth<-unique(RawGenotypes[,c("Cell","Pos","Depth")])
  Depthdic<-Depth$Depth
  names(Depthdic)<-paste(Depth$Cell, Depth$Pos,sep="")
  
  # Summarise
  Genotypes.summary<-table(paste(RawGenotypes$Cell,RawGenotypes$Variants,sep="_")) %>% as.data.frame()
  Genotypes.summary$Cell<-strsplit(as.character(Genotypes.summary$Var1),"_") %>% sapply(.,function(x){x[1]})
  Genotypes.summary$Variants<-strsplit(as.character(Genotypes.summary$Var1),"_") %>% sapply(.,function(x){paste(x[2:4],collapse="_")})
  Genotypes.summary$cellPos<-strsplit(as.character(Genotypes.summary$Var1),"_") %>% sapply(.,function(x){paste(x[1:2],collapse="")})
  Genotypes.summary$depth<-Depthdic[Genotypes.summary$cellPos]
  Genotypes.summary<-Genotypes.summary[,c("Var1","Cell","Variants","Freq","depth")]
  Genotypes.summary$Type<-strsplit(Genotypes.summary$Variants,"_") %>% sapply(.,function(x){paste(x[2],x[3],sep="_")})
  Genotypes.summary$Context<-ContextsDic[strsplit(Genotypes.summary$Variants,"_") %>% sapply(.,function(x){x[1]})]
  if(filterN){
    Genotypes.summary<-subset(Genotypes.summary,!grepl("N",Genotypes.summary$Context))
  }
  return(Genotypes.summary)
}

redeemR.read_edge_partition<-function(path,thr="S",Processed=F,rdsname="/VariantsGTSummary.RDS", remove_edge = TRUE, ratio_based_filter = TRUE){
  if(Processed){
    VariantsGTSummary<-readRDS(paste(path,"/VariantsGTSummary.RDS",sep=""))
  }else{
    if(missing(path)|missing(thr)){
      message("redeemR.read(path,thr)")
      message("missing variable path or thr")
      message("path is a string to the redeemV result folder that contains RawGenotypes.XX")
      message("thr is one of T,LS,S,VS:")
      message("T(Total),LS(Less Stringent:c=0.75,a1=2,a2=1), S(Stringent:c=0.75,a1=3,a2=2), VS(Very Stringent:c=0.75,a1=4,a2=3)")
      message("Term from redeemV is deprecated: VerySensitive equal to Less Stringent, Sensitive equal to Stringent, Specific equal to Very Stringent")
      return(NULL)
    }
    GiveName<-c("UMI","Cell","Pos","Variants","Call","Ref","FamSize","GT_Cts","CSS","DB_Cts","SG_Cts","Plus","Minus","Depth")
    if(thr=="T"){
      RawGenotypes<-read.table(paste(path,"/RawGenotypes.Total.StrandBalance",sep=""))
    }else if(thr=="LS"){
      RawGenotypes<-read.table(paste(path,"/RawGenotypes.VerySensitive.StrandBalance",sep=""))
    }else if(thr=="S"){
      RawGenotypes<-fread(paste(path,"/RawGenotypes.Sensitive.StrandBalance",sep="")) %>% data.frame()
    }else if(thr=="VS"){
      RawGenotypes<-read.table(paste(path,"/RawGenotypes.Specific.StrandBalance",sep=""))
    } 
    colnames(RawGenotypes)<-GiveName
    VariantsGTSummary<-GTSummary_CL_edge_partition(RawGenotypes, remove_edge = remove_edge, ratio_based_filter = ratio_based_filter)
    ##Calculate heteroplasmy
    VariantsGTSummary$hetero<-with(VariantsGTSummary,Freq/depth)
    attr(VariantsGTSummary,"thr")<-thr
    attr(VariantsGTSummary,"depth")<-DepthSummary(path)
    attr(VariantsGTSummary,"path")<-path
    saveRDS(VariantsGTSummary,paste(path,rdsname,sep=""))
    return(VariantsGTSummary)
  }
}

make_matrix_cl <- function(object,onlyhetero=T){
  require(dplyr)
  require(Matrix.utils)
  if(onlyhetero){
    GTsummary.filtered<-subset(object@GTsummary.filtered,!Variants %in% object@HomoVariants)
    message("Only heteroplasmic mutations are used")
  }
  Cts.Mtx<-dMcast(GTsummary.filtered,Cell~Variants,value.var = "Freq")
  colnames(Cts.Mtx)<-strsplit(as.character(colnames(Cts.Mtx)),"_") %>% sapply(.,function(x){paste(x[1],x[2],x[3],sep="")})
  
  # append cells that are missing
  missing_cells <- object@CellMeta$Cell[!(object@CellMeta$Cell %in% rownames(Cts.Mtx))]
  missing_mat <- sparseMatrix(
    i = 1:length(missing_cells),
    j = rep(dim(Cts.Mtx)[2], length(missing_cells)),
    x = 0
  )
  rownames(missing_mat) <- missing_cells; colnames(missing_mat) <- colnames(Cts.Mtx)
  Cts.Mtx_new <- rbind(Cts.Mtx, missing_mat)[object@CellMeta$Cell,]
  Cts.Mtx <- Cts.Mtx_new
  
  # back to normal
  Cts.Mtx.bi<-Cts.Mtx
  Cts.Mtx.bi[Cts.Mtx.bi>=1]<-1
  object@Cts.Mtx.bi<-Cts.Mtx.bi
  object@Cts.Mtx<-Cts.Mtx
  message("@Cts.Mtx and @Cts.Mtx.bi are added")
  #                 validObject(object)
  return(object)
}
