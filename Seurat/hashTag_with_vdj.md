# Read count data
**learn from "D:\study\singleCell\10Xtech\cihan_cell_hashing_slides_oct10_fnl.pdf"**
**"https://satijalab.org/seurat/articles/hashing_vignette.html"**


# hash tag
```
library(Seurat)
setwd("D:/work/project1/KRAS_control_VDJ_20211202/output")
mytcr_file = '../cellranger/KRAS_G12V_VDJ/outs/filtered_contig_annotations.csv'
myhto_file = '../cellranger/KRAS_G12V_hash/outs/raw_feature_bc_matrix'
n = 100
HTODemux_threshold = 0.99

# keep each cell sum(UMI) > n
mytcr <- read.csv(mytcr_file,header = T)# nrow(mytcr) [1] 25508
myhto <- Read10X(myhto_file) # ncol(myhto) 589964
myhto <- myhto[,colSums(myhto) > n] # ncol(myhto) [1] 22451

# remove cells not detected in scTCR-seq
keep_cells <- intersect(mytcr$barcode,colnames(myhto)) #length(keep_cells) #[1] 11115
myhto <- myhto[,keep_cells] #ncol(myhto) #[1] 11115

# create Seurat object
myobj <- CreateSeuratObject(counts = myhto, project = "G12V_hash") #dim(myobj)  2 11115
myfeature = as.data.frame(myobj@assays$RNA@counts)
head(myfeature[,1:4])
#                     AAACCTGAGCCAGTAG-1 AAACCTGAGCTGTTCA-1 AAACCTGCAAACCCAT-1 AAACCTGCAAACTGTC-1
#MHC-C0301-TotalSeqC                 12                 16                 38                162
#MHC-C0302-TotalSeqC                377                749                 50                124

### check boxplot count distribute
library(ggplot2)
library(reshape2)
pdata <- melt(as.matrix(myfeature))
ggplot(pdata,aes(x=Var1, y=value)) +
  geom_violin() +
  ylim(0,1000) + ylab("UMI") + xlab("")

# hash tag assignment
myobj@meta.data$mut = as.numeric(myfeature["MHC-C0301-TotalSeqC",])
myobj@meta.data$wt = as.numeric(myfeature["MHC-C0302-TotalSeqC",])

myobj[["HTO"]] <- CreateAssayObject(counts = myhto)
myobj <- NormalizeData(myobj, assay = "HTO", normalization.method = "CLR")
myobj <- HTODemux(myobj, assay = "HTO", kfunc = "clara", positive.quantile = HTODemux_threshold)

HTOHeatmap(myobj, assay = "HTO", ncells = 7600)
table(myobj@meta.data$hash.ID)
#            Doublet MHC-C0302-TotalSeqC MHC-C0301-TotalSeqC            Negative 
#                500                5234                4743                 638 

# check result by PLOT
# Group cells based on the max HTO signal
Idents(myobj) <- "HTO_maxID"
RidgePlot(myobj, assay = "HTO", features = rownames(myobj[["HTO"]])[1:2], ncol = 2)
Idents(myobj) <- "HTO_classification.global"
VlnPlot(myobj, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# output
write.csv(myobj@meta.data,"0_G12V_hashTag_result.csv",quote = F)
```

# scTCR-seq
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
setwd("D:/work/project1/KRAS_control_VDJ_20211202/output")
mytcr_file = '../cellranger/KRAS_G12V_VDJ/outs/filtered_contig_annotations.csv'
myhto_file = '1_G12V_hashTag_result.csv'
# read file 
mytcr <- read.csv(mytcr_file,header = T) #[1] 25508
myhto <- read.csv(myhto_file,row.names = 1) #[1] 11115
# TCR_Reshape
mytcr2 <- TCR_Reshape(mytcr)
# Tag_MergeTo_TCR
mytcr2 <- Tag_MergeTo_TCR(myhto,mytcr2)
# just keep hto useful cells
mytcr2$hash_sample <- mytcr2$hash.ID
mytcr2$hash_sample <- gsub("MHC-C0301-TotalSeqC","mut",mytcr2$hash_sample)
mytcr2$hash_sample <- gsub("MHC-C0302-TotalSeqC","wt",mytcr2$hash_sample)
mytcr2 <- mytcr2[mytcr2$hash_sample %in% c("mut","wt"),]
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
write.table(cloneFreq_abs,"1_G12V_myclonotype.txt",quote = F,row.names = F,sep = "\t")
write.table(mytcr3,"1_G12V_cell_TCR_meta.txt",quote = F,row.names = F,sep = "\t")
write.table(mydf,"1_G12V_cell_TCR_meta_long_format.txt",quote = F,row.names = F,sep = "\t")
```
# visualization
```
library(ggplot2)
library(reshape2)
setwd("D:/work/project1/KRAS_control_VDJ_20211202/output")
mydata <- read.delim("1_G12V_cell_TCR_meta.txt",header = T,row.names = 1,sep = "\t")
filter_freq = 30
mydata <- mydata[order(mydata$myclonotype_freq,decreasing = T),]

TCR_type <- unique(mydata[,c("myclonotype_abs","TRA_num","TRB_num")])
rownames(TCR_type) <- TCR_type$myclonotype_abs

mydata2 <- mydata[mydata$myclonotype_freq > filter_freq,]
length(table(mydata2$myclonotype_abs)) #[1] 31
pdata <- as.data.frame.matrix(table(mydata2$myclonotype_abs, mydata2$hash_sample))
pdata$mut_total = table(mydata$hash_sample)["mut"]
pdata$wt_total = table(mydata$hash_sample)["wt"]
pdata$mut_percentage <- round(pdata$mut / pdata$mut_total,3)
pdata$wt_percentage <- round(pdata$wt / pdata$wt_total,3)
pdata$TRATRB <- paste(paste0("(",TCR_type[rownames(pdata),]$TRA_num), 
                      paste0(TCR_type[rownames(pdata),]$TRB_num,")"),
                      sep = "+")

pval = c()
for(i in 1:nrow(pdata)){
  a = pdata[i,"mut"]
  b = pdata[i,"mut_total"] - a
  c = pdata[i,"wt"]
  d = pdata[i,"wt_total"] - c
  tmp <- matrix(c(a,b,c,d),nrow = 2)
  res = fisher.test(tmp)
  pval = c(pval,res$p.value)
}

pdata$fisher_pval = pval
pdata$fisher_sig = "ns"
pdata[pdata$fisher_pval >= 0.01 & pdata$fisher_pval < 0.05,]$fisher_sig = "*"
pdata[pdata$fisher_pval >= 0.001 & pdata$fisher_pval < 0.01,]$fisher_sig = "**"
pdata[pdata$fisher_pval < 0.001,]$fisher_sig = "***"

freq <- melt(as.matrix(pdata[,c("mut","wt")]),)
pdata_melt <- melt(as.matrix(pdata[,c("mut_percentage","wt_percentage")]),)
colnames(pdata_melt) <- c("myclonotype","sample","freq_percentage")
pdata_melt$freq = freq$value
pdata_melt$myclonotype <- factor(pdata_melt$myclonotype, levels = paste0("myclonotype",1:2000))

pdata_melt$freq = paste0("(",pdata_melt$freq,")")


p <- ggplot(pdata_melt,aes(x=myclonotype,y=freq_percentage,fill=sample)) +
  geom_bar(stat="identity", position=position_dodge(),width=0.7) +
  geom_text(aes(label = freq_percentage),
            #vjust=1.6,
            color="black",
            position = position_dodge(0.7), 
            size=3.5,
            hjust =  0,
            vjust = 0.5,
            angle = 90) +
  geom_text(aes(label = freq,y=freq_percentage+0.004),
            #vjust=1.6,
            color="black",
            position = position_dodge(0.7), 
            size=3.5,
            hjust =  0,
            vjust = 0.5,
            angle = 90) +
    scale_fill_manual(values=c('#E69F00','#999999')) +
  #facet_wrap(~sample, ncol=1) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.045)) +
  theme(
      #legend.position="none",
      text = element_text(size=15),
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 12),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color = "black")) +
  ggtitle("G12V_TCR_clone_frequency > 30") +
  xlab("") +
  ylab("Frequency") +
  coord_cartesian(clip='off')

for (i in 1:nrow(pdata)) {
  p <- p + annotation_custom(
    textGrob(
      label=pdata[paste0("myclonotype",i),]$TRATRB, 
      rot=45, gp=gpar(fontsize=10,col="blue")),
      xmin= i+0.1, xmax= i+0.1, ymin=-0.005, ymax=0
    )
}


for (i in 1:nrow(pdata)) {
  p <- p + annotation_custom(
    textGrob(
      label=pdata[paste0("myclonotype",i),]$fisher_sig, 
      rot=45, gp=gpar(fontsize=10,col="red")),
      xmin= i-0.6, xmax= i-0.6, ymin=-0.01, ymax=0
    )
}

pdf("2_G12V_PLOT_TCR_clone_frequency_more_30.pdf",height = 7,width = 14)
p
dev.off()

### PLOT plan2
pdata_melt$TRATRB <- pdata[pdata_melt$myclonotype,]$TRATRB
pdata_melt$fisher_sig <- pdata[pdata_melt$myclonotype,]$fisher_sig
pdata_melt$info <- paste(pdata_melt$myclonotype, paste(pdata_melt$fisher_sig,pdata_melt$TRATRB,sep = "---"))
addline_format <- function(x,...){
    gsub('\\s','\n',x)
}

ggplot(pdata_melt,aes(x=myclonotype,y=freq_percentage,fill=sample)) +
  geom_bar(stat="identity", position=position_dodge(),width=0.7) +
  geom_text(aes(label = freq_percentage),
            #vjust=1.6,
            color="black",
            position = position_dodge(0.7), 
            size=3.5,
            hjust =  0,
            vjust = 0.5,
            angle = 90) +
  geom_text(aes(label = freq,y=freq_percentage+0.004),
            #vjust=1.6,
            color="black",
            position = position_dodge(0.7), 
            size=3.5,
            hjust =  0,
            vjust = 0.5,
            angle = 90) +
    scale_fill_manual(values=c('#E69F00','#999999')) +
  #facet_wrap(~sample, ncol=1) +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.045)) +
  scale_x_discrete(labels = addline_format(pdata_melt$info)) +
  theme(
      #legend.position="none",
      text = element_text(size=15),
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 12),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color = "black")) +
  ggtitle("G12V_TCR_clone_frequency > 30") +
  xlab("") +
  ylab("Frequency") 
```

# G12V freq > 10 & fisher.test(tmp,alternative = “greater”) sig = “***”  
clonotype frequency >= 10  
Mut-percentage >> WT-percentage  
Just TRA-TRB or TRA1-TRA2-TRB clonotype
```
library(ggplot2)
library(reshape2)
setwd("D:/work/project1/KRAS_control_VDJ_20211202/output")
mydata <- read.delim("1_G12V_cell_TCR_meta.txt",header = T,row.names = 1,sep = "\t")
filter_freq = 10
mydata <- mydata[order(mydata$myclonotype_freq,decreasing = T),]

TCR_type <- unique(mydata[,c("myclonotype_abs","TRA_num","TRB_num")])
rownames(TCR_type) <- TCR_type$myclonotype_abs

mydata2 <- mydata[mydata$myclonotype_freq >= filter_freq,]
length(table(mydata2$myclonotype_abs)) #[1] 31
pdata <- as.data.frame.matrix(table(mydata2$myclonotype_abs, mydata2$hash_sample))
pdata$mut_total = table(mydata$hash_sample)["mut"]
pdata$wt_total = table(mydata$hash_sample)["wt"]
pdata$mut_percentage <- round(pdata$mut / pdata$mut_total,3)
pdata$wt_percentage <- round(pdata$wt / pdata$wt_total,3)
pdata$TRATRB <- paste(paste0("(",TCR_type[rownames(pdata),]$TRA_num), 
                      paste0(TCR_type[rownames(pdata),]$TRB_num,")"),
                      sep = "+")

pval = c()
for(i in 1:nrow(pdata)){
  a = pdata[i,"mut"]
  b = pdata[i,"mut_total"] - a
  c = pdata[i,"wt"]
  d = pdata[i,"wt_total"] - c
  tmp <- matrix(c(a,b,c,d),nrow = 2)
  res = fisher.test(tmp,alternative = "greater")
  pval = c(pval,res$p.value)
}

pdata$fisher_pval = pval
pdata$fisher_sig = "ns"
pdata[pdata$fisher_pval >= 0.01 & pdata$fisher_pval < 0.05,]$fisher_sig = "*"
pdata[pdata$fisher_pval >= 0.001 & pdata$fisher_pval < 0.01,]$fisher_sig = "**"
pdata[pdata$fisher_pval < 0.001,]$fisher_sig = "***"

vals <- as.numeric(gsub("myclonotype","",rownames(pdata)))
pdata <- pdata[order(vals),]
pdata <- pdata[pdata$fisher_sig == "***",]
```
# extract sequencing function
```
library(ggplot2)
library(reshape2)
library(plyr)
library(Biostrings)
setwd('D:/work/project1/KRAS_control_VDJ_20211202/output')
# Input file
input_meta = '1_G12V_cell_TCR_meta_long_format.txt'
input_fasta <- "../cellranger/KRAS_G12V_VDJ/outs/consensus.fasta"

my_clonotypes <- rownames(pdata)
mymeta <- read.table(input_meta,header = T,sep = "\t")
topclone <- mymeta[mymeta$myclonotype_abs %in% my_clonotypes &
                     mymeta$TRA_num >= 1 & mymeta$TRB_num == 1,
                     c("chain","v_gene","j_gene","cdr3",
                       "raw_clonotype_id","raw_consensus_id","myclonotype_abs",
                       "myclonotype_freq","TRA_num","TRB_num","chain2")]
dim(topclone)
#[1] 2532    6

mydata <- topclone[!duplicated(topclone[,c("raw_consensus_id","myclonotype_abs")]),]
dim(mydata)
#[1] 105   6

fdta <- readDNAStringSet(input_fasta)
df <- data.frame(fdta)
colnames(df) = "Sequence"
df$consensus_ID <- names(fdta)
rownames(df) <- names(fdta)
df <- df[rownames(df) %in% mydata$raw_consensus_id,,drop=F] #dim(df) 105

setdiff(mydata$raw_consensus_id,rownames(df))

mydata$sequence <- df[mydata$raw_consensus_id,]$Sequence
mydata <- mydata[!duplicated(mydata[,c("chain","cdr3","myclonotype_abs","chain2")]),]
mydata <- mydata[,c(7:10,1:6,11,12)]
vals <- as.numeric(gsub("myclonotype","", mydata$myclonotype_abs))
mydata <- mydata[order(vals),]

myfinal <- merge(mydata,pdata,
                 by.x = "myclonotype_abs", by.y = 0, all.x = T)
myfinal <- myfinal[,c(1:10,13:21,11,12)]
myfinal$mutPercent_dividedBy_wtPercent <- myfinal$mut_percentage / myfinal$wt_percentage
myfinal <- myfinal[,c(1:16,22,17:21)]
write.table(myfinal, file = "3_G12V_selected_clonotype_sequence.txt",row.names = F,quote = F,sep = "\t")
# clonotype frequency >= 10 and Mut-percentage >> WT-percentage 
```
