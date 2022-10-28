library(TCGAbiolinks)
library(SummarizedExperiment)
#BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data") 
#BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

# Genome of reference: hg38

##----------------- download clinical data -----------------
project <- 'TCGA-LAML'
type <- 'Clinical'
#download_dir <- 'D:/work/task_other/AML/TCGA_LAML'

clinicals <- GDCquery_clinic(
  project = project,
  type = type
  )
  
##----------------- download SNV data -----------------

project <- 'TCGA-LAML'
data_category <- 'Simple Nucleotide Variation'
data_type <- 'Masked Somatic Mutation'
access <- 'open'
legacy <- FALSE
#download_dir <- 'D:/work/task_other/AML/TCGA_LAML'

query <- GDCquery(
  project = project, 
  data.category = data_category, 
  access = access, 
  legacy = legacy, 
  data.type = data_type
  )

GDCdownload(query)
maf <- GDCprepare(query)

##----------------- download expression data -----------------
project <- 'TCGA-LAML'
data_category <- 'Transcriptome Profiling'
data_type <- 'Gene Expression Quantification'
workflow_type <- 'STAR - Counts'
#download_dir <- 'D:/work/task_other/AML/TCGA_LAML'

query <- GDCquery(
  project = project,
  data.category = data_category,
  data.type = data_type, 
  workflow.type = workflow_type
  )
GDCdownload(query)
exp <- GDCprepare(query)

# extract counts and tmp
counts <- assay(exp,i = "unstranded") #tpm_unstrand fpkm_unstrand
tpm <- assay(exp,i = "tpm_unstrand")
geneid22=data.frame(id=rowData(exp)@listData[["gene_id"]], 
                    gene_name= rowData(exp)@listData[["gene_name"]],
                    gene_type=rowData(exp)@listData[["gene_type"]])

counts=cbind(geneid22,counts)
tpm=cbind(geneid22,tpm)
#a1=gsub("_PAR_Y","",counts$id)
#a2=gsub("_PAR_Y","",tpm$id)

##----------------- save raw data from GDCquery -----------------
saveRDS(list(clinicals = clinicals, maf = maf, counts = counts, tpm = tpm),
        "D:/work/task_other/AML/TCGA_LAML/TCGA_LAML_Clinical_Maf_Exp.rds")


