# use Biostrings package to read or write fasta files

```R
library(ggplot2)
library(Biostrings)
# library(seqinr)
library(data.table)
main_path = 'D:/work/'
setwd(main_path)
clone_meta_file = 'D:/work/VDJ/consensus_annotations.csv'
myVDJ_ref_file = 'D:/work/reference/cellranger/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0/fasta/regions.fa'

# defined clone type 
verify_clone_pos <- c(2,3,5,6,9,7,11,4,17,31,8,35)
verify_clone_neg <- c(1,12,13,131)
verify_df <- data.frame(clonotype = paste0('clonotype', 
                                           c(verify_clone_pos,verify_clone_neg)),
                        verify = c(rep('pos', length(verify_clone_pos)),
                                   rep('neg', length(verify_clone_neg))))
rownames(verify_df) <- verify_df$clonotype
#                 clonotype verify
# clonotype2     clonotype2    pos
# clonotype3     clonotype3    pos
# clonotype5     clonotype5    pos

# extract selected clonotype Vgenes
clone_meta <- read.csv(clone_meta_file)
clone_meta <- clone_meta[clone_meta$clonotype_id %in% verify_df$clonotype,]
clone_meta$verify <- verify_df[clone_meta$clonotype_id,]$verify

#----------------------------------------------
#          read referene fasta file 
#----------------------------------------------
ref <- readDNAStringSet(myVDJ_ref_file)
head(ref)
# DNAStringSet object of length 6:
#     width seq                                                                   names               
# [1]  1031 AGTCTGCGAGAAATCCCACCATCTACCCACTGA...GTGATCATGTCAGAGGGAGATGGCATCTGCTAC 1|IGHA ENSMUST000...
# [2]  1166 AGTCTGCGAGAAATCCCACCATCTACCCACTGA...GGCCCGTTTGGCAGCAAAGAGGTCCCCCAGTAC 2|IGHA ENSMUST000...
# [3]   773 GTAATGAAAAGGGACCTGACATGTTCCTCCTCT...CCTTCCAGGAGACCTGATGGTCCTGCCCTTGCC 3|IGHD ENSMUST000...
# [4]   872 GTAATGAAAAGGGACCTGACATGTTCCTCCTCT...TACAGTGGCTTCGTCACCTTCATCAAGGTGAAG 4|IGHD ENSMUST000...

ref <- as.data.frame(ref)
ref <- ref[!(rownames(ref) %like% "5'UTR"),,drop = F]
ref <- ref[!(rownames(ref) %like% "GRCm38-release94"),,drop = F]
ref <- ref[!(rownames(ref) %like% "IGH"),,drop = F]
names <- rownames(ref)
tmp <- strsplit(names,'\\|')
tmp <- as.data.frame(do.call(rbind,tmp))
ref$names <- rownames(ref)
rownames(ref) <- tmp$V3

pos_fa <- ref[unique(clone_meta[clone_meta$verify == "pos",]$v_gene),,drop=F]
neg_fa <- ref[unique(clone_meta[clone_meta$verify == "neg",]$v_gene),,drop=F]

pos_fa <- pos_fa[order(rownames(pos_fa)),,drop=F]
neg_fa <- neg_fa[order(rownames(neg_fa)),,drop=F]

#----------------------------------------------
#          data.frame to Biostrings format
#----------------------------------------------
pos_fa2 = DNAStringSet(pos_fa$x)
names(pos_fa2) <- pos_fa$names
neg_fa2 = DNAStringSet(neg_fa$x)
names(neg_fa2) <- neg_fa$names

# clone_meta$ID <- ref[clone_meta$v_gene,]$names
# write.csv(clone_meta, "my_scTCRseq_verify_pos_neg_meta.csv",
#          row.names = F,quote = F)
          
#----------------------------------------------
#          write fasta/fastq file to file
#----------------------------------------------
writeXStringSet(pos_fa2, "tumorSpecificVerify_positive_clonotype_Vgenes.fa", 
                format="fasta")
writeXStringSet(neg_fa2, "tumorSpecificVerify_negtive_clonotype_Vgenes.fa", 
                format="fasta")

```


