library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(circlize)
workpath <- 'D:/work/task_other/AML/TCGA_LAML/'
data_file <- 'TCGA_LAML_Clinical_Maf_Exp.rds'
MutGene <- c('NPM1','CEBPA','KIT','DNMT3A','IDH1','IDH2','NARS1','TET2','ASXL1','RUNX1','TP53')
ReceptorGene <- c("MICA","MICB","ULBP1","ULBP2","ULBP3","RAET1E","RAET1G","RAET1L")

# read Clinical, MAf, Count and TPM data
setwd(workpath)
datalist <- readRDS(data_file)

#-----------------------------------------------
# maf reshape
#-----------------------------------------------
maf <- datalist$maf
a <- str_split(maf$Tumor_Sample_Barcode,'-',simplify = TRUE)
maf$Sample <- paste(a[,1], a[,2], a[,3], sep = "-")
maf <- maf[,c("Hugo_Symbol","Entrez_Gene_Id","Variant_Classification",
              "Variant_Type","dbSNP_RS","Sample")]

#-----------------------------------------------
# tpm reshape
#-----------------------------------------------

tpm <- datalist$tpm
tpm <- tpm[tpm$gene_name %in% ReceptorGene,]
rownames(tpm) <- tpm$gene_name
tpm <- tpm %>% select(-c(1:3))
a <- str_split(colnames(tpm),'-',simplify = TRUE)
tpm_samples <- paste(a[,1], a[,2], a[,3], sep = "-")
colnames(tpm) <- tpm_samples

#-----------------------------------------------
# Just keep exp and maf both true samples 
#-----------------------------------------------
keep_samples <- intersect(colnames(tpm),maf$Sample) # 99 patients
maf <- maf[maf$Sample %in% keep_samples,]
tpm <- tpm[,keep_samples]

#-------------------------------------------------
# Maf PLOT data prapare 
#-------------------------------------------------
# Maf file Keep selected MutGene
mymaf <- maf[maf$Hugo_Symbol %in% MutGene,
             c("Hugo_Symbol","Variant_Classification","Sample")]
mymaf <- unique(mymaf)

# Maf reshape to wide
mymaf_wide <- mymaf %>%
  pivot_wider(names_from = Sample, values_from = Variant_Classification)
mymaf_wide <- as.data.frame(mymaf_wide)
rownames(mymaf_wide) <- mymaf_wide$Hugo_Symbol
mymaf_wide <- mymaf_wide %>% select(-Hugo_Symbol)

# convert each mymaf_wide[i,j] to character sperate by ";"
for (i in 1:nrow(mymaf_wide)) {
  for (j in 1:ncol(mymaf_wide)) {
    mymaf_wide[i,j] = paste(mymaf_wide[i,j][[1]], collapse = ";")
  }
}

# add non-Mutation patient samples
NonMutSample <- setdiff(keep_samples,colnames(mymaf_wide))
tmp <- data.frame(matrix(ncol = length(NonMutSample), nrow = nrow(mymaf_wide)))
colnames(tmp) <- NonMutSample
mymaf_wide <- cbind(mymaf_wide,tmp)
mymaf_wide[is.na(mymaf_wide)] <- ""

# convert to matrix
mymaf_wide <- as.matrix(mymaf_wide)

table(mymaf$Variant_Classification)
# Frame_Shift_Del   Frame_Shift_Ins      In_Frame_Del      In_Frame_Ins Missense_Mutation Nonsense_Mutation       Splice_Site 
#               4                13                 1                 3                26                 2                 4

#----------------------------------------------------------
#  PLOT
#----------------------------------------------------------

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#d9d9d9", col = NA))
  },
  # bug red
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),  
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"),h*0.33, 
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["In_Frame_Ins"], col = NA))
  },
  # small green
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  #
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33,  
              gp = gpar(fill = col["Splice_Site"], col = NA))
  }
)

col = c("Frame_Shift_Del" = "#f0a7bf", "Frame_Shift_Ins" = "#1e2e75", 
        "In_Frame_Del" = "#f38d22","In_Frame_Ins" = "#624698", 
        "Missense_Mutation" = "#e72520", "Nonsense_Mutation" = "#19572b",
        "Splice_Site" = "#64b040")

heatmap_legend_param = list(title = "Alternations", 
                            at = c("Missense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                                   "In_Frame_Del","In_Frame_Ins",
                                   "Nonsense_Mutation","Splice_Site"))
#pdf("1.pdf",onefile =F,height = 5,width = 9)
p <- oncoPrint(mymaf_wide,
          alter_fun = alter_fun, 
          col = col, 
          #row_order = 1:nrow(wide), column_order = sample_order,
          #column_title = "Altered in 133(96%) of 138 samples",
          heatmap_legend_param = heatmap_legend_param
)
#dev.off()


#----------------------------------------------------------
#   PLOT add expression heatmap: order by MutationCondition
#----------------------------------------------------------
oncoPrint_order1 <- colnames(mymaf_wide)[unname(p@column_order)]
LogTpm1 <- LogTpm[,oncoPrint_order1]
mymaf_wide1 <- mymaf_wide[,oncoPrint_order1]
mycolor2 = colorRamp2(c(0,2,3,4,5,6), c("#f7fbff","#c6dbef","#9ecae1","#6baed6","#2171b5","#08306b"))

pdf("1_order_by_MutationCondition.pdf",height = 5, width = 15)
oncoPrint(mymaf_wide1,
          alter_fun = alter_fun, 
          col = col, 
          column_order = oncoPrint_order1,
          heatmap_legend_param = heatmap_legend_param) %v%
  Heatmap(as.matrix(LogTpm1), 
          cluster_columns = F,
          column_names_gp = gpar(fontsize = 8),
          row_names_gp = gpar(fontsize = 10),
          name = "Log2(TPM+1)", 
          col = mycolor2,
          height = unit(3, "cm"))
dev.off()

#----------------------------------------------------------
#   PLOT add expression heatmap: order by expression heatmap 
#----------------------------------------------------------
p_order <- pheatmap(LogTpm[c("MICB","MICA"),],clustering_method = "average")
oncoPrint_order2 <- p_order$tree_col$labels[p_order$tree_col$order]
LogTpm2 <- LogTpm[,oncoPrint_order2]
mymaf_wide2 <- mymaf_wide[,oncoPrint_order2]
mycolor2 = colorRamp2(c(0,2,3,4,5,6), c("#f7fbff","#c6dbef","#9ecae1","#6baed6","#2171b5","#08306b"))

pdf("2_order_by_expressionLevels.pdf",height = 5, width = 15)
oncoPrint(mymaf_wide2,
          alter_fun = alter_fun, 
          col = col, 
          column_order = oncoPrint_order2,
          heatmap_legend_param = heatmap_legend_param) %v%
  Heatmap(as.matrix(LogTpm2), 
          cluster_columns = F,
          column_names_gp = gpar(fontsize = 8),
          row_names_gp = gpar(fontsize = 10),
          name = "Log2(TPM+1)", 
          col = mycolor2,
          height = unit(3, "cm"))
dev.off()

