
## ref
Protocol  
https://www.biolegend.com/en-us/totalseq/single-cell-rna  
https://assets.ctfassets.net/an68im79xiti/6OUafQzYFi6cPH8pqKzYF/2e498570823168843887da238a3f86b0/CG000391_CellLabelingwithCellMultiplexingOligo_RevA.pdf

Analysis  
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis

10x genome VDJ sequencing  
http://www.sci666.net/51441.html  
https://cloud.tencent.com/developer/article/1635484

https://cloud.tencent.com/developer/article/1814314  #velocyto  
https://cloud.tencent.com/developer/article/1819555  #Garnett

Total seq  
http://www.bio-city.net/index.php/NewsProducts/conts/id/4136  
https://www.biolegend.com/en-us/search-results/totalseq-c0301-anti-mouse-hashtag-1-antibody-17157?GroupID=GROUP20  
https://www.biolegend.com/en-us/products/totalseq-c0302-anti-mouse-hashtag-2-antibody-17158?GroupID=GROUP20

feature barcode  
What is Feature Barcode Technology for Immune Profiling?  
https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/feature-bc  
scRNA-seq https://adinasarapu.github.io/posts/2019/01/blog-post-sc-ranseq/

-------------- plan 1 -----------------------------------  
step1: Feature barcode only analysis - cellranger count  
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/no-gex-analysis  
To use cellranger in Feature Barcode Only mode, follow instructions for Feature Barcode Analysis, and omit Gene Expression entries from the Libraries CSV file.  
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis  
step2: cellranger vdj  
normal vdj analysis methods  
step3: barcode mapping

-------------- plan 2 -----------------------------------  
Analyzing V(D)J and Gene Expression / Feature Barcode with cellranger multi  
https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/multi  
https://support.10xgenomics.com/single-cell-vdj/software/analysis-of-multiple-libraries/latest/overview  
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi  
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_cp#csv

```shell
cellranger multi --id=Jurkat_Raji_10K --csv=SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K.csv
```
### SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K.csv
```
[gene-expression]
ref,/path/to/refdata-gex-GRCh38-2020-A
expect-cells,10000

[libraries]
fastq_id,fastqs,lanes,feature_types
SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_1_gex,/path/to/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_1_gex,any,Gene Expression
SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_1_multiplexing_capture,/path/to/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_1_multiplexing_capture,any,Multiplexing Capture
SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_2_multiplexing_capture,/path/to/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_2_multiplexing_capture,any,Multiplexing Capture

[samples]
sample_id,cmo_ids,description
Jurkat,CMO301,Jurkat
Raji,CMO302,Raji
```
### 3' Gene Expression with Feature Barcoding and Cell Multiplexing
```
[gene-expression]
reference,/path/to/transcriptome

[libraries]
fastq_id,fastqs,feature_types
gex1,/path/to/fastqs,Gene Expression
abc1,/path/to/fastqs,Antibody Capture
mux1,/path/to/fastqs,Multiplexing Capture

[samples]
sample_id,cmo_ids
sample1,CMO301|CMO302
sample2,CMO303|CMO304
```

