#!/bin/bash

dir0=/mnt/data1/John/Pioneer_Factors

cd $dir0/Analysis/ATAC_seq/Union/R_output_Common_Merge0/Clusters
for i in `ls *15-clusters.bed`; do
  filename=${i%-clusters.bed}
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' $i > ${filename}_5_column.bed ##Adds a unique peak name and strand info
  cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed .
  mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
  bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b ${filename}_5_column.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > ${filename}_peaks_subtract.bed
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${filename}_peaks_subtract.bed > ${filename}_peaks_subtract_5_column.bed ##Adds a unique peak name and strand info
  findMotifsGenome.pl ${filename}_5_column.bed mm10 ${filename}_motif/ -size 200 -p 35 -bg ${filename}_peaks_subtract_5_column.bed -keepFiles	##Does motif analysis at unique ATAC seq peaks
done

cd $dir0/Analysis/ATAC_seq/Union/R_output_Common_1.5/Clusters
for i in `ls *clusters.bed`; do
  filename=${i%-clusters.bed}
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' $i > ${filename}_5_column.bed ##Adds a unique peak name and strand info
  cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed .
  mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
  bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b ${filename}_5_column.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > ${filename}_peaks_subtract.bed
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${filename}_peaks_subtract.bed > ${filename}_peaks_subtract_5_column.bed ##Adds a unique peak name and strand info
  #findMotifsGenome.pl ${filename}_5_column.bed mm10 ${filename}_motif/ -size 200 -p 35 -bg ${filename}_peaks_subtract_5_column.bed -keepFiles   ##Does motif analysis at unique ATAC seq peaks
done
