#!/bin/bash
#Homer motif analysis for common CD4, CD8, and NK peaks


cd ${dir0}/Analysis/Motif/Homer/Common_CD4_CD8_NK
#---------------------CREATE 6 COLUMN BED FILE--------------------------------------------
##This generates a 6 column bedgraph that includes unique peak names and adds strand info used for homer
awk 'BEGIN {OFS="\t"} {print $1, $2, $3}' $dir0/Analysis/ATAC_seq/R_output/CD4_CD8_NK_common_B_Lsk_FC.bedgraph > CD4_CD8_NK_common_B_Lsk_FC_5_column_temp.bed ##Adds a unique peak name and strand info

awk 'BEGIN {OFS="\t"} {print $0 "\tpeak" NR "\t0" "\t+"}' CD4_CD8_NK_common_B_Lsk_FC_5_column_temp.bed > CD4_CD8_NK_common_B_Lsk_FC_5_column.bed ##Adds a unique peak name and strand info
batchFindMotifsGenome.pl mm10 -size 1000 -p 35 -f CD4_CD8_NK_common_B_Lsk_FC_5_column.bed	##Does motif analysis at unique ATAC seq peaks