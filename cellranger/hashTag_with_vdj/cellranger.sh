## Just HTO + TCR
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


## HTO + RNA + TCR
#------------------- colitisBrainColon -------------------------------
# PLAN 1 ## cell ranger nulti
~/tools/cellranger-6.1.2/cellranger multi --id=cBrain_HTO_VDJ_GEX \
	--csv=../script/cellranger_multi_config_cBrain.csv \
	--localcores=8 \
	--localmem=32
	
# PLAN 2 ## sample hash
~/tools/cellranger-6.1.2/cellranger count --id=CBrain_hash \
                                          --libraries=../script/cellranger_meta/library_CBrain.csv \
                                          --transcriptome=/home/innovent/myproject/reference/10x/refdata-gex-mm10-2020-A \
                                          --feature-ref=../script/cellranger_meta/featureRef_CBrain.csv
#### VDJ
~/tools/cellranger-6.1.2/cellranger vdj --id=CBrain_VDJ \
                                        --reference=/home/innovent/myproject/reference/10x/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0 \
                                        --fastqs=/home/innovent/myproject/Colitis_20211206/rawdata/colitis_SPLN \
                                        --sample=colitisBrainColon_T_20211126NA_TTGCCCGTGC \
                                        --localcores=8 \
                                        --localmem=32



