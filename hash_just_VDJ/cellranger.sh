#------------------- G12V -------------------------------
#### sample hash
~/tools/cellranger-6.1.2/cellranger count --id=KRAS_G12V_hash \
                                          --libraries=../script/cellranger_meta/library_G12V.csv \
                                          --transcriptome=/home/innovent/myproject/reference/10x/refdata-gex-mm10-2020-A \
                                          --feature-ref=../script/cellranger_meta/featureRef_G12V.csv
#### VDJ
~/tools/cellranger-6.1.2/cellranger vdj --id=KRAS_G12V_VDJ \
                                        --reference=/home/innovent/myproject/reference/10x/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0 \
                                        --fastqs=/home/innovent/myproject/KRAS_control_VDJ_20211202/rawdata/G12V \
                                        --sample=A2_G12V_T_20211112NA_CACAATCCCA \
                                        --localcores=8 \
                                        --localmem=32


#------------------- G12D -------------------------------
#### sample hash
~/tools/cellranger-6.1.2/cellranger count --id=KRAS_G12D_hash \
                                          --libraries=../script/cellranger_meta/library_G12D.csv \
                                          --transcriptome=/home/innovent/myproject/reference/10x/refdata-gex-mm10-2020-A \
                                          --feature-ref=../script/cellranger_meta/featureRef_G12D.csv
#### VDJ
~/tools/cellranger-6.1.2/cellranger vdj --id=KRAS_G12D_VDJ \
                                        --reference=/home/innovent/myproject/reference/10x/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0 \
                                        --fastqs=/home/innovent/myproject/KRAS_control_VDJ_20211202/rawdata/G12D \
                                        --sample=A2_G12D_T_20211112NA_TCCGGGACAA \
                                        --localcores=8 \
                                        --localmem=32


#------------------- G12V # 20211116 hash-------------------------------
#### sample hash 
##### library_G12V_20211116.csv 
~/tools/cellranger-6.1.2/cellranger count --id=KRAS_G12V_hash_20211116 \
                                          --libraries=../script/cellranger_meta/library_G12V_20211116.csv \
                                          --transcriptome=/home/innovent/myproject/reference/10x/refdata-gex-mm10-2020-A \
                                          --feature-ref=../script/cellranger_meta/featureRef_G12V.csv



#-------------------- cellranger multi method --------------------------
### output is difficult to follow
~/tools/cellranger-6.1.2/cellranger multi --id=KRAS_G12V_hash_muli \
       	--csv=/home/innovent/myproject/KRAS_control_VDJ_20211202/script/cellranger_meta/cellranger_multi_config_G12V.csv
~/tools/cellranger-6.1.2/cellranger multi --id=KRAS_control_G12D \
	--csv=/home/innovent/myproject/KRAS_control_VDJ_20211202/script/cellranger_meta/cellranger_multi_config_G12D.csv


~/tools/cellranger-6.1.2/cellranger multi --id=KRAS_control_G12V_20211116 \
        --csv=/home/innovent/myproject/KRAS_control_VDJ_20211202/script/cellranger_meta/cellranger_multi_config_G12V_20211116.csv
