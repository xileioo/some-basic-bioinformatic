# igBlast_Alpaca_library

https://ncbi.github.io/igblast/cook/How-to-set-up.html <br/>

IgBLAST internally only supports organisms including human, mouse, rat, rabbit and rhesus_monkey. However, starting with version 1.16.0, IgBLAST supports custom organism, although you need to prepare additional files as follows. Additionally, starting with version 1.20.0 of igblastn, you can also use the custom FWR/CDR annotation function for custom organism (see procedure to use custom FWR/CDR annotation below).

Make the germline V gene annotation file for your own organism. The purpose of this file is to provide information such as CDR/FWR regions for your germline V genes (see internal_data/human/human.ndm.imgt for an example). This file must be put under the internal_data/my_organism folder (for example internal_data/sheep). For annotations using the IMGT numbering system, the file must be named my_organism.ndm.imgt for igblastn (for example sheep.ndm.imgt) or my_organism.pdm.imgt for igblastp (for example sheep.pdm.imgt). If you use the Kabat numbering system, the file must be named my_organism.pdm.kabat or my_organism.pdm.kabat. Note that you do not need to make annotations for all alleles of a V gene. You only need to annotate one allele.

Make a blast sequence database for germline V genes that correspond to what you have annotated above. This database needs to be named my_organism_V (for example sheep_V). Make sure you use the -parse_seqids flag when using makeblastdb. The blast database files need to be put under internal_data/my_organism folder (for example internal_data/sheep). Note that this database is intended only as internal data for IgBLAST and does not need to be updated unless there is a new germline V gene (not new allele). Typically, the germline V gene database you want to search (i.e., specified by -germline_db_V parameter) is a different one (for example, the one that contains all alleles).

If you also want CDR3/FWR4 information, you need to supply a file that has information such as CDR3 stop for germline J genes (see optional_file/human_gl.aux for an example). This file can have any name and can be put anywhere as long as you supply it to the -auxiliary_data parameter when running IgBLAST(for example -auxiliary_data my_foler/my_file).

To run IgBLAST for your organism, please make sure you specify the -organism my_organism parameter.

List of documents to be prepared：<br/>
1. alpaca.ndm.imgt:FWR1+CDR1+FWE2+CDR2+FWR3
2. blast sequence database: Makeblastdb
3. auxiliary_data: CDR3+FWR4

Vicugna pacos download site: https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Vicugna_pacos/IG/ <br/>
```
$ perl ~/tools/igblast/igblast.1.21.0/bin/edit_imgt_file.pl IGHV.fasta > IGHV_2.fasta
$ perl ~/tools/igblast/igblast.1.21.0/bin/edit_imgt_file.pl IGHD.fasta > IGHD_2.fasta
$ perl ~/tools/igblast/igblast.1.21.0/bin/edit_imgt_file.pl IGHJ.fasta > IGHJ_2.fasta
$ makeblastdb -parse_seqids -dbtype nucl -in IGHV_2.fasta -out blastn_refdb/alpaca_V
$ makeblastdb -parse_seqids -dbtype nucl -in IGHD_2.fasta -out blastn_refdb/alpaca_D
$ makeblastdb -parse_seqids -dbtype nucl -in IGHJ_2.fasta -out blastn_refdb/alpaca_J
```

Commond <br/>
```
~/tools/igblast/igblast.1.21.0/bin/igblastn \n
-germline_db_V ~/tools/igblast/IMGT_IG_reference_alpaca/blastn_refdb/alpaca_V \
-germline_db_D ~/tools/igblast/IMGT_IG_reference_alpaca/blastn_refdb/alpaca_D \
-germline_db_J ~/tools/igblast/IMGT_IG_reference_alpaca/blastn_refdb/alpaca_J \
-organism alpaca \
-domain_system imgt  \
-auxiliary_data ~/tools/igblast/igblast.1.21.0/optional_file/alpaca_gl.aux \
-ig_seqtype Ig \
-show_translation \
-outfmt 19 \
-num_threads 8 \
-query flash_output/*.extendedFrags.fasta \
-out igblast_output/*.igblast.fmt19.out > log/*.igblast.log 2>&1
```
