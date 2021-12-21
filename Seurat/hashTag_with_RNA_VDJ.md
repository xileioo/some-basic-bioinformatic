# Seurat object and UMAP
```
library(Seurat)
setwd("D:/work/output")
mytcr_file = '../cellranger/HTO_VDJ_GEX/outs/per_sample_outs/HTO_VDJ_GEX/vdj_t/filtered_contig_annotations.csv'
myhto_file = '../cellranger/HTO_VDJ_GEX/outs/per_sample_outs/HTO_VDJ_GEX/count/sample_feature_bc_matrix'
Count_HTO_threshold = 80
HTODemux_threshold = 0.98

# createAssay
mycount <- Read10X(myhto_file)
myobj  = CreateSeuratObject(counts = mycount$`Gene Expression`,project = "cBrain",min.cells = 3)
myobj[["HTO"]] <- CreateAssayObject(counts = mycount$`Antibody Capture`) 

# filter nCount_HTO < nCount_HTO_threshold
myobj[["percent.mt"]] <- PercentageFeatureSet(myobj, assay = "RNA",pattern = "^mt-")
VlnPlot(myobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
myobj <- subset(myobj, subset = nCount_HTO > Count_HTO_threshold)
myobj <- subset(myobj, subset = nFeature_RNA > 200 & nCount_RNA < 15000 & percent.mt < 5)
dim(myobj)#[1] 15735 13008

#  remove cells not detected in scTCR-seq
mytcr <- read.csv(mytcr_file,header = T)
keep_cells <- intersect(mytcr$barcode,colnames(myobj)) # length(keep_cells) #[1] 9600
myobj <- subset(myobj, cells = keep_cells)
dim(myobj)
# [1] 15735  9600
VlnPlot(myobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","nCount_HTO"), ncol = 4)
VlnPlot(myobj, features = c("nCount_HTO"), ncol = 1) + ylim(1,500)

# hash tag assignment
myobj <- NormalizeData(myobj, assay = "HTO", normalization.method = "CLR")
myobj <- HTODemux(myobj, assay = "HTO", kfunc = "clara", positive.quantile = HTODemux_threshold)

HTOHeatmap(myobj, assay = "HTO", ncells = 7600)
table(myobj@meta.data$hash.ID)
#Doublet            Negative MHC-C0302-TotalSeqC MHC-C0301-TotalSeqC MHC-C0304-TotalSeqC MHC-C0303-TotalSeqC 
#   2762                3086                2426                 325                 718                 283 

# check result by PLOT
# Group cells based on the max HTO signal
Idents(myobj) <- "HTO_maxID"
RidgePlot(myobj, assay = "HTO", features = rownames(myobj[["HTO"]])[1:4], ncol = 2)
Idents(myobj) <- "HTO_classification.global"
VlnPlot(myobj, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# myobj@meta.data$con_Brain = as.numeric(myobj@assays$HTO@data["MHC-C0301-TotalSeqC",])
# myobj@meta.data$con_Colon = as.numeric(myobj@assays$HTO@data["MHC-C0302-TotalSeqC",])
# myobj@meta.data$ATB_Brain = as.numeric(myobj@assays$HTO@data["MHC-C0303-TotalSeqC",])
# myobj@meta.data$ATB_Colon = as.numeric(myobj@assays$HTO@data["MHC-C0304-TotalSeqC",])
myobj@meta.data$hash_sample <- myobj@meta.data$hash.ID
myobj@meta.data$hash_sample <- gsub("MHC-C0301-TotalSeqC","Con_Brain",myobj@meta.data$hash_sample)
myobj@meta.data$hash_sample <- gsub("MHC-C0302-TotalSeqC","Con_Colon",myobj@meta.data$hash_sample)
myobj@meta.data$hash_sample <- gsub("MHC-C0303-TotalSeqC","ATB_Brain",myobj@meta.data$hash_sample)
myobj@meta.data$hash_sample <- gsub("MHC-C0304-TotalSeqC","ATB_Colon",myobj@meta.data$hash_sample)
# output
write.csv(myobj@meta.data,"0_cBrain_hashTag_result.csv",quote = F)

#-------------------------------------------------------------------------------
# Cell cluster
#------------------------------------------------------------------------------- 
myobj <- subset(myobj,subset = hash_sample %in% c("Con_Brain","Con_Colon","ATB_Brain","ATB_Colon"))
myobj <- NormalizeData(myobj, normalization.method = "LogNormalize", scale.factor = 10000)
myobj <- FindVariableFeatures(myobj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(myobj), 10)
#[1] "Gzma"      "Il17a"     "Il22"      "Ifitm1"    "Cxcl3"     "Penk"      "Hspa1a"    "Hist1h2ae" "Hspa1b"
myobj <- ScaleData(myobj, verbose = FALSE)
myobj <- RunPCA(myobj, npcs = 30, verbose = FALSE)
myobj <- RunUMAP(myobj, reduction = "pca", dims = 1:30)
myobj <- FindNeighbors(myobj, reduction = "pca", dims = 1:30)
myobj <- FindClusters(myobj, resolution = 0.5)
DefaultAssay(myobj) <- "RNA"
#myobj <- FindClusters(myobj, resolution = 0.3)
# Visualization
DimPlot(myobj, reduction = "umap")
DimPlot(myobj, reduction = "umap",group.by = "hash_sample")
FeaturePlot(myobj,features = c("Cd4","Cd8a","Pdcd1","Tigit"))
FeaturePlot(myobj,features = c("Gzmk","Gzmb","Itgae","Ifngr1"))
#-------------------------------------------------------------------------------
# Cell type assignment
#------------------------------------------------------------------------------- 
saveRDS(myobj, file = "1_cell_UMAP.rds")
```

# TCR clonotype defined
```
TCR_Reshape <- function(mytcr){
  mytcr <- mytcr[mytcr$cdr3 != "None",]
  mytcr2 <- data.frame(barcode = character(),
                     chain2 = character(),
                     chain_nt = character(),
                     raw_clonotype_id = character())
  
  for(barcode in unique(mytcr$barcode)){
    #barcode = "AAACCTGCAAACCCAT-1"
    mydf <- mytcr[mytcr$barcode == barcode,-c(12:22,25,26)]
    mydf <- mydf[order(mydf$chain,mydf$v_gene,mydf$j_gene,mydf$cdr3,decreasing = F),]
    chain = ""
    chain_nt = ""
    for(j in 1:nrow(mydf)){
      info = paste(mydf[j,]$v_gene,mydf[j,]$j_gene,mydf[j,]$cdr3,sep = ",")
      each_chain = paste(mydf[j,]$chain,info,sep = ":")
      each_chain_nt = paste(mydf[j,]$chain,mydf[j,]$cdr3_nt,sep = ":")
      chain = paste(chain,each_chain,sep = ";")
      chain_nt = paste(chain_nt,each_chain_nt,sep = ";")
      }
    chain = substring(chain, 2) # remove first sep (;) from previous step
    chain_nt = substring(chain_nt, 2)
    raw_clonotype_id = unique(mydf$raw_clonotype_id)
    tmp <- data.frame(barcode = barcode,
                      chain2 = chain,
                      chain_nt = chain_nt,
                      raw_clonotype_id = raw_clonotype_id)
    mytcr2 <- rbind(mytcr2,tmp)
  }
  rownames(mytcr2) <- mytcr2$barcode
  return(mytcr2)
}

Tag_MergeTo_TCR <- function(myhto,mytcr2){
  mytcr2 <- merge(mytcr2,myhto[,c("nCount_RNA","hash.ID")],all.x = T,
             by.x = "barcode", by.y = 0)
  #res_htoNum <- table(mytcr2[!duplicated(mytcr2$barcode),]$hash.ID)
  return(mytcr2)
}

TCR_MycloneDefined_Hto <- function(mytcr2){
  cloneFreq_abs <- as.data.frame(table(mytcr2$chain2))
  cloneFreq_abs <- cloneFreq_abs[order(cloneFreq_abs$Freq, decreasing = T),]
  cloneFreq_abs$myclonotype_abs <- paste("myclonotype",1:nrow(cloneFreq_abs),sep = "")
  colnames(cloneFreq_abs) = c("chain2","myclonotype_freq","myclonotype_abs")
  cloneFreq_abs$chain2 <- as.character(cloneFreq_abs$chain2)
  return(cloneFreq_abs)
}

Dual_TCR_Defined <- function(cloneFreq_abs){
  cloneFreq_abs$TRA_num = 0
  cloneFreq_abs$TRB_num = 0
  for (i in 1:nrow(cloneFreq_abs)){
    TRA_num <- length(grep("TRA",strsplit(cloneFreq_abs[i,]$chain2,";")[[1]],value = T))
    TRB_num <- length(grep("TRB",strsplit(cloneFreq_abs[i,]$chain2,";")[[1]],value = T))
    cloneFreq_abs[i,]$TRA_num = TRA_num
    cloneFreq_abs[i,]$TRB_num = TRB_num
  }
  return(cloneFreq_abs)
}


### 
library(dplyr)
setwd("D:/work/project1/Colitis_20211206/output")
mytcr_file = '../cellranger/HTO_VDJ_GEX/outs/per_sample_outs/HTO_VDJ_GEX/vdj_t/filtered_contig_annotations.csv'
myhto_file = '0_cBrain_hashTag_result.csv'
# read file 
mytcr <- read.csv(mytcr_file,header = T) #nrow(mytcr) [1] 28778
myhto <- read.csv(myhto_file,row.names = 1) #nrow(myhto) [1] 9600
# TCR_Reshape
mytcr2 <- TCR_Reshape(mytcr) # dim(mytcr2) [1] 15092     4
# Tag_MergeTo_TCR
mytcr2 <- Tag_MergeTo_TCR(myhto,mytcr2)
# just keep hto useful cells
mytcr2$hash_sample <- mytcr2$hash.ID
mytcr2$hash_sample <- gsub("MHC-C0301-TotalSeqC","Con_Brain",mytcr2$hash_sample)
mytcr2$hash_sample <- gsub("MHC-C0302-TotalSeqC","Con_Colon",mytcr2$hash_sample)
mytcr2$hash_sample <- gsub("MHC-C0303-TotalSeqC","ATB_Brain",mytcr2$hash_sample)
mytcr2$hash_sample <- gsub("MHC-C0304-TotalSeqC","ATB_Colon",mytcr2$hash_sample)

mytcr2 <- mytcr2[mytcr2$hash_sample %in% c("Con_Brain","Con_Colon","ATB_Brain","ATB_Colon"),]
dim(mytcr2)
# [1] 3752    7

# Myclonotype defined
cloneFreq_abs <- TCR_MycloneDefined_Hto(mytcr2)
cloneFreq_abs <- Dual_TCR_Defined(cloneFreq_abs)
# merge mytcr2 and cloneFreq_abs
mytcr3 <- merge(mytcr2,cloneFreq_abs,by="chain2",all.x = T)
# merge mytc2 and mytcr
mydf <- merge(mytcr, 
              mytcr3[,c("barcode","hash.ID","hash_sample","myclonotype_abs","myclonotype_freq",
                                            "TRA_num","TRB_num","chain2")],
              by = "barcode")
# Output
mytcr3 <- mytcr3[,c("barcode","hash.ID","hash_sample","nCount_RNA","raw_clonotype_id",
                    "myclonotype_abs","myclonotype_freq","TRA_num","TRB_num",
                    "chain2","chain_nt")]
write.table(cloneFreq_abs,"1_cBrain_myclonotype.txt",quote = F,row.names = F,sep = "\t")
write.table(mytcr3,"1_cBrain_cell_TCR_meta.txt",quote = F,row.names = F,sep = "\t")
write.table(mydf,"1_cBrain_cell_TCR_meta_long_format.txt",quote = F,row.names = F,sep = "\t")
```

