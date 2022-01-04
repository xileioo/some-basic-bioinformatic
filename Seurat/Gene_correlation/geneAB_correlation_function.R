# correlation between geneA and geneB
# calulation from Seurat object of single cell RNA-seq data 
# by lixue 20210903 

# ---------------------- example --------------------------------------
# library(Seurat)
# source("geneAB_correlation_function.R") 
# input_file = '1_CD8pos_cells_normalization_notlog_GSE144945.rds'
# myobj <- readRDS(input_file) # read seurat normalized data
# myobj <- subset(x = myobj, subset = PDCD1 > 0 & GZMB > 0) # keep cells with geneA>0 & geneB>0
# GENE_COR_PLOT is a function to calculate the correlation between geneA and geneB
# GENE_COR_PLOT(myobj = myobj, 
#               group = "patient", # group is the colum name from seurat@meta.data 
#               geneA = "PDCD1",  # gene symbol - A
#               geneB = "GZMB", # gene symbol -B
#               datasetID="GSE144945", # Any string
#               tumorType ="NSCLC") # Any string
# ggsave(filename = "1_cor_PD1_GZMB_GSE144945_NSCLC.pdf",height = 8, width = 16)
# ggsave(filename = "1_cor_PD1_IL2RA_GSE144945_NSCLC.png",height = 6, width = 12,device = "png")
# ---------------------- example --------------------------------------



# library
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(EnvStats)


# !(x %in% y)
'%!in%' <- function(x,y)!('%in%'(x,y))

# GET_METADATA function: get metadata and the expresion level(RNA assay) of geneA and geneB
GET_METADATA <- function(myobj, geneA, geneB) {
  exp <- GetAssayData(object = myobj, slot = 'data',assay="RNA")
  if (geneA %!in% rownames(exp)) {
    print(paste("There is no ",geneA,sep = ""))
    return(NA)
  }else if (geneB %!in% rownames(exp)) {
    print(paste("There is no ",geneB,sep = ""))
    return(NA)
  }
  mymeta <- myobj@meta.data
  mymeta$geneA <- exp[geneA,]
  mymeta$geneB <- exp[geneB,]
  return(mymeta)
}

# PLOT_SCATTER function: scatter plot, geneA is x-axis, geneB is y-axis
PLOT_SCATTER <- function(mymeta, group, datasetID="test", tumorType ="test",geneA, geneB){
  # arguments
  plot_title = paste(datasetID,tumorType,sep = "-") #plot title
  text_fontSize = 5 
  axis_fontSize = 15
  max_x_Decile = max(mymeta$geneA)/10
  max_y_Decile = max(mymeta$geneB)/10
  mycolors = c(brewer.pal(9,"Paired"),brewer.pal(8,"Dark2"))
  cells_anno <- paste("total cell = ", nrow(mymeta), sep = "")
  # plot
  p1 <- ggplot(data=mymeta, aes(x=geneA, y=geneB)) + 
    geom_point(data=mymeta, aes(x=geneA, y=geneB, color=get(group)),size = 2,alpha = 0.8) +
    scale_color_manual(values = mycolors) +
    stat_cor(method = 'pearson',size = text_fontSize,
             label.x = max_x_Decile*6, 
             label.y = max_y_Decile*8) +
    stat_cor(method = 'spearman',size = text_fontSize,
             label.x = max_x_Decile*6, 
             label.y = max_y_Decile*7) +
    annotate(geom="text", label="pearson", size=text_fontSize,
             x=max_x_Decile*5, y=max_y_Decile*8) +
    annotate(geom="text", label="spearman",size=text_fontSize,
             x=max_x_Decile*5, y=max_y_Decile*7) +
    annotate(geom="text", label=cells_anno,size=text_fontSize,
             x=max_x_Decile*6, y=max_y_Decile*9) +
    theme_bw() +
    ggtitle(plot_title) +
    xlab(paste(geneA,"expression level",sep = " ")) + 
    ylab(paste(geneB,"expression level",sep = " ")) +
    theme(panel.grid = element_blank(),
        #legend.position = c(),
        legend.title = element_text(size = axis_fontSize),
        axis.text=element_text(size=axis_fontSize, colour = "black"),
        axis.title=element_text(size=axis_fontSize),
        plot.title = element_text(size=axis_fontSize),
        legend.text = element_text(size=axis_fontSize))  +
    guides(color = guide_legend(override.aes = list(size = 5)))
  return(p1)
}

# PLOT_VIOLIN function: violin plot, geneA is group into high and low expression levels
PLOT_VIOLIN <- function(mymeta, datasetID="test", tumorType ="test",geneA, geneB){
  # arguments
  plot_title = paste(datasetID,tumorType,sep = "-")
  text_fontSize = 5
  axis_fontSize = 15
  # plot
  mymeta$geneA_level <- paste(geneA,"_high",sep = "")
  mymeta[mymeta$geneA <= median(mymeta$geneA),]$geneA_level <- paste(geneA,"_low",sep = "")
  p2 <- ggplot(data = mymeta, aes(x=geneA_level, y=geneB, color=geneA_level)) +
    geom_violin(size = 2) +
    geom_jitter(size = 1,alpha = .5) +
    scale_color_manual(values = c("#357FBC","#45AD4C")) +
    stat_compare_means(method = "wilcox.test",
                     method.args = list(alternative = "greater"),
                     size = text_fontSize) +
    stat_n_text(size = text_fontSize) +
    theme_bw() +
    xlab(paste(geneA,"expression level",sep = " ")) + 
    ylab(paste(geneB,"expression level",sep = " ")) +
    theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text=element_text(size=axis_fontSize, colour = "black"),
        axis.title=element_text(size=axis_fontSize),
        plot.title = element_text(size=axis_fontSize)) +
    ggtitle(plot_title) 
  return(p2)
}

# pipeline od correlation plot
GENE_COR_PLOT <- function(myobj,group="patient",geneA,geneB,datasetID="test", tumorType ="test"){
  mymeta <- GET_METADATA(myobj = myobj, geneA = geneA, geneB = geneB)
  P1 <- PLOT_SCATTER(mymeta,group = group, datasetID, tumorType,geneA = geneA, geneB = geneB)
  P2 <- PLOT_VIOLIN(mymeta,datasetID,tumorType,geneA = geneA, geneB = geneB)
  P3 = P1 + P2 # merge P1 and P2 into P3
  return(P3)
}
