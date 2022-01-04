### simple merge count single cell RNA-seq

```R
library(Seurat)
# Working directory
main_path = 'd:/work/project2/PMID_34290408/downloadData/GSE176021/mRNA/tumor/'
setwd(main_path)
# Input path
mrna_tumor_path <- "d:/work/project2/PMID_34290408/downloadData/GSE176021/mRNA/tumor/"
tcr_tumor_path <- "d:/work/project2/PMID_34290408/downloadData/GSE176021/VDJ/table/"
patient_meta <- "D:/work/project2/PMID_34290408/downloadData/patient_meta.csv"# Output path
out_path = 'D:/work/project2/PMID_34290408/output/'


## get total tumor matrix file folder
mrna_folder <- list.files(mrna_tumor_path)
tcr_files <- list.files(tcr_tumor_path)
## just seleted "Squamous Cell Carcinoma" samples
patient_meta <- read.csv(patient_meta,header = T,row.names = 1)
scc_patient_ID <- rownames(patient_meta[patient_meta$Histology == "Squamous Cell Carcinoma",])
mrna_folder <- mrna_folder[grep(paste(scc_patient_ID,collapse = "|"), mrna_folder)]


i = 1  
for (each_sample in mrna_folder) {
  #each_sample = "MD01-005_tumor_2"
  myobj <- Read10X(data.dir = each_sample)
  #myobj <- myobj[,myobj["CD4",]==0 & myobj["CD8A",]>0]
  myobj <- myobj[,myobj["CD8A",] > myobj["CD4",]]
  #myobj <- myobj[,myobj["PDCD1",]>0]
  myobj <- CreateSeuratObject(counts = myobj, project = each_sample, min.cells = 3, min.features = 200)
  myobj[["percent.mt"]] <- PercentageFeatureSet(myobj, pattern = "^MT-")
  myobj[["percent.ribo"]] <- PercentageFeatureSet(myobj,
                                                  pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  #VlnPlot(myobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 2)
  myobj <- subset(myobj, subset = nFeature_RNA > 200 & percent.mt < 10 & percent.ribo > 10)
  
  ### filter no TCR cells
  mytcr_file <- grep(each_sample, tcr_files, value = T)
  if(identical(mytcr_file, character(0))){
    print(paste("there is no TCR files for sample",each_sample))
    next}
  mytcr <- read.csv(file = paste("../../VDJ/table/",mytcr_file,sep = ""),header = T,row.names = 1)
  cells_use <- intersect(rownames(myobj@meta.data), rownames(mytcr))
  myobj <- subset(myobj, cells = cells_use)
  myobj@meta.data <- cbind(myobj@meta.data, mytcr[rownames(myobj@meta.data),])
  if (i == 1) {
    myobj.combined <- myobj
  } else {
    myobj.combined <- merge(myobj.combined, y = myobj)
  }
  i = i + 1
  print(i)
}
#dim(myobj.combined)
saveRDS(myobj.combined, file = paste(out_path,"1_total_Tcells_CD8morethanCD4_scc_mergeCount.rds",sep = ""))
```
### simple merge single cell RNA-seq - UMAP
```R
library(Seurat)
# Working directory
main_path = 'D:/work/project2/PMID_34290408/output/'
setwd(main_path)
# Input path
input_data <- "1_total_Tcells_CD8morethanCD4_scc_mergeCount.rds"# Output path
output_prefix = '1_total_Tcells_CD8morethanCD4_scc_'


myobj <- readRDS(input_data)
dim(myobj) #[1] 15896 30964


tmp <- strsplit(as.character(myobj@meta.data$orig.ident),'_tumor')
tmp <- do.call(rbind, tmp)
myobj@meta.data$patient <- tmp[,1]


myobj <- NormalizeData(myobj, normalization.method = "LogNormalize", scale.factor = 10000)
myobj <- FindVariableFeatures(myobj, selection.method = "vst", nfeatures = 2000)
myobj <- ScaleData(myobj, features = rownames(myobj))
myobj <- RunPCA(myobj, npcs = 20, verbose = FALSE)
myobj <- RunUMAP(myobj, reduction = "pca", dims = 1:20)
myobj <- FindNeighbors(myobj, reduction = "pca", dims = 1:20)
myobj <- FindClusters(myobj, resolution = 0.5)
#myobj <- FindClusters(myobj, resolution = 0.3)# Visualization
saveRDS(myobj,paste(output_prefix,"simpleMerge_UMAP.rds",sep = ""))

DimPlot(myobj, reduction = "umap")
DimPlot(myobj, reduction = "umap",group.by = "patient")
FeaturePlot(myobj,features = c("CD4","CD8A","PDCD1","GZMB"))
```
### Integration single cell RNA-seq - UMAP
```R
library(Seurat)
# Working directory
main_path = 'D:/work/project2/PMID_34290408/output/'
setwd(main_path)
# Input path
input_data <- "1_total_Tcells_CD8morethanCD4_scc_mergeCount.rds"# Output path
output_prefix = '1_total_Tcells_CD8morethanCD4_scc_'


myobj <- readRDS(input_data)
dim(myobj) #[1] 15896 30964


#myobj@meta.data$patient <- substr(myobj@meta.data$orig.ident,1,nchar(myobj@meta.data$orig.ident)-2)
tmp <- strsplit(as.character(myobj@meta.data$orig.ident),'_tumor')
tmp <- do.call(rbind, tmp)
myobj@meta.data$patient <- tmp[,1]


# split the dataset into a list
myobj.list <- SplitObject(myobj, split.by = "patient")
myobj.list <- lapply(X = myobj.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = myobj.list, nfeatures = 2000)
myobj.list <- PrepSCTIntegration(object.list = myobj.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = myobj.list,
                                         normalization.method = "SCT",
                                         anchor.features = features)


memory.limit()
memory.limit(size=10^9)
myobj <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT") ## memory limit
myobj <- RunPCA(myobj, verbose = FALSE)
myobj <- RunUMAP(myobj, reduction = "pca", dims = 1:20)
myobj <- FindNeighbors(myobj, reduction = "pca", dims = 1:20)
myobj <- FindClusters(myobj, resolution = 0.5)
DefaultAssay(myobj) <- "SCT"
dim(myobj)
#[1] 13447  8709


saveRDS(myobj,paste(output_prefix,"Integrate_UMAP.rds",sep = ""))


DimPlot(myobj, reduction = "umap", group.by = "patient")
DimPlot(myobj, reduction = "umap", label = TRUE, repel = TRUE)
FeaturePlot(myobj,features = c("CD4","CD8A","PDCD1","CD3E"))
FeaturePlot(myobj,features = c("SELL","CCR7","CD44","IL7R"))
```

# script 2
```R
# 20210022 # Seurat3 - MD01-005,MD01-004,MD043-011

library(Seurat)
library(data.table)
# Working directory
main_path = 'D:/work/project2/PMID_34290408/downloadData/GSE176021/mRNA/tumor'
setwd(main_path)
# Input
scrna_dir = c("MD01-004_tumor_1","MD01-005_tumor_2","MD01-005_tumor_3",
              "MD01-005_tumor_4","MD01-005_tumor_5","MD01-005_tumor_6",
              "MD01-005_tumor_7","MD01-005_tumor_8","MD01-005_tumor_9",
              "MD043-011_tumor_2","MD043-011_tumor_3","MD043-011_tumor_4",
              "MD043-011_tumor_5")
#Output path
out_path = 'd:/work/project2/PMID_34290408/output'


#---------------------------------------------------------------------
# loop create Seurat Object
#---------------------------------------------------------------------
for (i in 1:length(scrna_dir)) {
  scrna_dir_one <- scrna_dir[i]
  mycount <- Read10X(data.dir = scrna_dir_one)
  myobj <- CreateSeuratObject(counts = mycount, project = scrna_dir_one, 
                              min.cells = 3, min.features = 250) 
  myobj[["percent.mt"]] <- PercentageFeatureSet(myobj, pattern = "^MT-")
  # myobj[["percent.ribo"]] <- PercentageFeatureSet(myobj, 
  #                                               pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  VlnPlot(myobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  myobj <- subset(myobj, subset = percent.mt < 10)
  
  genes <- rownames(myobj)
  genes <- genes[-(which(genes %like% "^MT-"))] # remove MT genes
  genes <- genes[-(which(genes %like% "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA"))] # remove ribo genes
  genes <- genes[-(which(genes %like% "^TRA|^TRB|^TRD|^TRG"))] # remove TRA/TRB/TRD/TRG genes
  myobj <- subset(myobj, features = genes)
  myobj@meta.data <- cbind(samples = rownames(myobj@meta.data),myobj@meta.data)
  if (i == 1) {
    myobj.merge <- myobj
  }else{
    myobj.merge <- merge(myobj.merge,myobj)
  }
  print(scrna_dir_one)
  print(dim(myobj))
}

dim(myobj.merge) # [1] 18320 99812
```

#### simply merge
```R
### simple merge single cell RNA-seq - UMAP
myobj <- myobj.merge
tmp <- strsplit(as.character(myobj@meta.data$orig.ident),'_tumor')
tmp <- do.call(rbind, tmp)
myobj@meta.data$patient <- tmp[,1]

myobj <- NormalizeData(myobj, normalization.method = "LogNormalize", scale.factor = 10000)
myobj <- FindVariableFeatures(myobj, selection.method = "vst", nfeatures = 3000)
myobj <- ScaleData(myobj, features = rownames(myobj))
myobj <- RunPCA(myobj, npcs = 20, verbose = FALSE)
myobj <- RunUMAP(myobj, reduction = "pca", dims = 1:20)
myobj <- FindNeighbors(myobj, reduction = "pca", dims = 1:20)
myobj <- FindClusters(myobj, resolution = 0.5)
#myobj <- FindClusters(myobj, resolution = 0.3)# Visualization
saveRDS(myobj,paste(output_prefix,"simpleMerge_UMAP.rds",sep = ""))

DimPlot(myobj, reduction = "umap")
DimPlot(myobj, reduction = "umap",group.by = "patient")
FeaturePlot(myobj,features = c("CD4","CD8A","PDCD1","GZMB"))
```
