## Read count data
**learn from "D:\study\singleCell\10Xtech\cihan_cell_hashing_slides_oct10_fnl.pdf"**
**"https://satijalab.org/seurat/articles/hashing_vignette.html"**
```R
library(Seurat)
setwd("D:/work/project1/KRAS_G12V_G12D/output/")
mycount <- Read10X("D:/work/project1/KRAS_G12V_G12D/KRAS/outs/filtered_feature_bc_matrix")
myobj <- CreateSeuratObject(counts = mycount, project = "KRAS_featureBarcode")
myfeature = as.data.frame(myobj@assays$RNA@counts)
head(myfeature[,1:4])
#                     AAACCTGAGATGCGAC-1 AAACCTGAGGAACTGC-1 AAACCTGAGTTGAGAT-1
# MHC-C0301-TotalSeqC                 19                293                574
# MHC-C0302-TotalSeqC               1454                143                103
myobj@meta.data$G12D = as.numeric(myfeature["MHC-C0301-TotalSeqC",])
myobj@meta.data$G12V = as.numeric(myfeature["MHC-C0302-TotalSeqC",])

myobj[["HTO"]] <- CreateAssayObject(counts = mycount)
myobj <- NormalizeData(myobj, assay = "HTO", normalization.method = "CLR")
myobj <- HTODemux(myobj, assay = "HTO", kfunc = "clara", positive.quantile = 0.99)

HTOHeatmap(myobj, assay = "HTO", ncells = 8000)
table(myobj@meta.data$hash.ID)
  # Doublet MHC-C0302-TotalSeqC MHC-C0301-TotalSeqC 
  #     669                3897                2893 
write.csv(myobj@meta.data,"1_KRAS_hashTag_result.csv",quote = F)
```
## TCR data annotation
```R
library(dplyr)
setwd("D:/work/project1/KRAS_G12V_G12D/output/")
mytcr <- read.csv("../KRAS_VDJ/outs/filtered_contig_annotations.csv",header = T)# 17484 
mytag <- read.csv("1_KRAS_hashTag_result.csv",row.names = 1)

cells = unique(intersect(rownames(mytag), mytcr$barcode))
length(cells)
# [1] 5906

mytcr$tag = "Nontag"
mytcr[mytcr$barcode %in% cells,]$tag <- mytag[mytcr[mytcr$barcode %in% cells,]$barcode,]$hash.ID

# > length(unique(mytcr[mytcr$tag == "MHC-C0301-TotalSeqC",]$barcode))
# [1] 2463
# > length(unique(mytcr[mytcr$tag == "MHC-C0302-TotalSeqC",]$barcode))
# [1] 2982
# > length(unique(mytcr[mytcr$tag == "Non_tag",]$barcode))
# [1] 2565
# > length(unique(mytcr[mytcr$tag == "Doublet",]$barcode))
# [1] 461

#-------------------------------------------------------------------------------
# TCR data reshape
#------------------------------------------------------------------------------- 
mytcr <- mytcr[mytcr$cdr3 != "None",] #[1] 17484 
mytcr2 <- data.frame(barcode = character(),
                     chain = character(),
                     chain_nt = character(),
                     raw_clonotype_id = character())

for(i in unique(mytcr$barcode)){
  #i = 'AAAGTAGTCTGGTTCC-1'
  mydf <- mytcr[mytcr$barcode == i,-c(12:22,25,26)]
  mydf <- mydf[order(mydf$chain,mydf$cdr3,decreasing = F),]
  chain = ""
  chain_nt = ""
  for(j in 1:nrow(mydf)){
    each_chain = paste(mydf[j,]$chain,mydf[j,]$cdr3,sep = ":")
    each_chain_nt = paste(mydf[j,]$chain,mydf[j,]$cdr3_nt,sep = ":")
    chain = paste(chain,each_chain,sep = ";")
    chain_nt = paste(chain_nt,each_chain_nt,sep = ";")
  }
  chain = substring(chain, 2) # remove first sep (;) from previous step
  chain_nt = substring(chain_nt, 2)
  barcode = i
  raw_clonotype_id = unique(mydf$raw_clonotype_id)
  tmp <- data.frame(barcode = barcode,
                    chain = chain,
                    chain_nt = chain_nt,
                    raw_clonotype_id = raw_clonotype_id)
  mytcr2 <- rbind(mytcr2,tmp)
}

dim(mytcr2) # [1] 8471    4
rownames(mytcr2) <- mytcr2$barcode

#-------------------------------------------------------------------------------
# Cell tag assignment
#------------------------------------------------------------------------------- 
mytcr2$tag = "NonTag"
mytcr2[mytcr2$barcode %in% cells,]$tag <- mytag[mytcr2[mytcr2$barcode %in% cells,]$barcode,]$hash.ID
#-------------------------------------------------------------------------------
# clone type assignment
#-------------------------------------------------------------------------------
# clone_freq base filter cells
clone_freq <- as.data.frame(table(mytcr2$raw_clonotype_id))
clone_freq <- clone_freq[order(clone_freq$Freq, decreasing = T),]
rownames(clone_freq) <- clone_freq$Var1
# > head(clone_freq)
#            Var1 Freq
# 1    clonotype1   51
# 1112 clonotype2   36
# 2223 clonotype3   31
# 3334 clonotype4   30
# 4445 clonotype5   19
# 5556 clonotype6   17
mytcr2$raw_cloneFreq <- clone_freq[mytcr2$raw_clonotype_id,]$Freq

# cloneFreq_absolute
cloneFreq_abs <- mytcr2[,c("chain","tag")] %>% 
                        group_by(tag,chain) %>% 
                        summarise(Freq=n())
cloneFreq_abs <- as.data.frame(cloneFreq_abs)
cloneFreq_abs <- cloneFreq_abs[order(cloneFreq_abs$tag, 
                                     cloneFreq_abs$Freq, decreasing = T),]
cloneFreq_abs$myclonotype_abs = "non"
total_num = table(cloneFreq_abs$tag)
cloneFreq_abs[cloneFreq_abs$tag == "NonTag",]$myclonotype_abs = paste("myclonotype",1:total_num["NonTag"],sep = "")
cloneFreq_abs[cloneFreq_abs$tag == "Doublet",]$myclonotype_abs = paste("myclonotype",1:total_num["Doublet"],sep = "")
cloneFreq_abs[cloneFreq_abs$tag == "MHC-C0301-TotalSeqC",]$myclonotype_abs = paste("myclonotype",1:total_num["MHC-C0301-TotalSeqC"],sep = "")
cloneFreq_abs[cloneFreq_abs$tag == "MHC-C0302-TotalSeqC",]$myclonotype_abs = paste("myclonotype",1:total_num["MHC-C0302-TotalSeqC"],sep = "")

# cloneFreq_absolute_name
cloneFreq_abs$sample <- cloneFreq_abs$tag
cloneFreq_abs[cloneFreq_abs$tag == "MHC-C0301-TotalSeqC",]$sample = "G12D"
cloneFreq_abs[cloneFreq_abs$tag == "MHC-C0302-TotalSeqC",]$sample = "G12V"
table(cloneFreq_abs$sample)
# Doublet    G12D    G12V  NonTag 
#     457    2293    2458    2292 
cloneFreq_abs$myclonotype_abs2 = paste(cloneFreq_abs$sample,
                                       cloneFreq_abs$myclonotype_abs,sep = "_")

TCR_share <- as.data.frame.matrix(table(cloneFreq_abs$chain,cloneFreq_abs$sample)) 
TCR_share$chain = rownames(TCR_share)
cloneFreq_abs <- merge(cloneFreq_abs,TCR_share,by=c("chain"))
cloneFreq_abs$share_TCR <- "ignore"
cloneFreq_abs[cloneFreq_abs$G12D == 1 & cloneFreq_abs$G12V == 1,]$share_TCR <- "TCR_sharein_G12V_G12D" 
# share clone is low frequency clone - meaning nothing

cloneFreq_abs <- cloneFreq_abs[order(cloneFreq_abs$sample, 
                                     cloneFreq_abs$Freq, decreasing = T),]

#-------------------------------------------------------------------------------
# dual TCR assignment
#------------------------------------------------------------------------------- 
cloneFreq_abs$TRA_num = 0
cloneFreq_abs$TRB_num = 0
for (i in 1:nrow(cloneFreq_abs)) {
  TRA_num <- length(grep("TRA",strsplit(cloneFreq_abs[i,]$chain,";")[[1]],value = T))
  TRB_num <- length(grep("TRB",strsplit(cloneFreq_abs[i,]$chain,";")[[1]],value = T))
  cloneFreq_abs[i,]$TRA_num = TRA_num
  cloneFreq_abs[i,]$TRB_num = TRB_num
}

table(cloneFreq_abs$TRA_num)
#   0    1    2 
# 885 5620  995
table(cloneFreq_abs$TRB_num)
#   0    1    2 
# 191 6689  620 

cloneFreq_abs$TRA_num <- as.character(cloneFreq_abs$TRA_num)
cloneFreq_abs[cloneFreq_abs$TRA_num == '1',]$TRA_num = "single"
cloneFreq_abs[cloneFreq_abs$TRA_num == '2',]$TRA_num = "dual"
table(cloneFreq_abs$TRA_num)
#   0   dual single 
# 885    995   5620
cloneFreq_abs$TRB_num <- as.character(cloneFreq_abs$TRB_num)
cloneFreq_abs[cloneFreq_abs$TRB_num == '1',]$TRB_num = "single"
cloneFreq_abs[cloneFreq_abs$TRB_num == '2',]$TRB_num = "dual"
table(cloneFreq_abs$TRB_num)
#   0   dual single 
# 191    620   6689 

#-------------------------------------------------------------------------------
# save file
#------------------------------------------------------------------------------- 

write.csv(cloneFreq_abs,"2_KRAS_TCR_cloneFreq_meta.csv",quote = F,row.names = F)

mytcr3 <- merge(mytcr2,cloneFreq_abs,by=c("tag","chain"))
mytcr3 <- mytcr3[,c(3,9,1,5,6,8,7,10:17,2,4)]
rownames(mytcr3) = mytcr3$barcode
mytcr3 <- mytcr3[,-1]
write.csv(mytcr3,"3_KRAS_cell_TCR_meta.csv",quote = F)
```

