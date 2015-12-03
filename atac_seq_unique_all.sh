#!/bin/bash

#This script will take all the union bedgraph file for all the ATAC-seq peaks generated earlier (atac_union_coverage.sh)
#Run them in an R script that will generate UNIQUE Peaks for T and NK cells
#Do motif analysis on these regions
cd $dir0/Analysis/ATAC_seq/Union/R_output

#---------------------R Script that gives ATAC UNIQUE PEAKS-------------------------------
#Rscript --verbose ${r_dir}/Unique_ATAC_seq_Peaks.R $dir0	

#-----------Copy Union Peak File Log 2----------------------------------------------------
cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed $dir0/Analysis/ATAC_seq/Union/R_output_FC1/
#cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed $dir0/Analysis/ATAC_seq/Union/R_output_Common/Intersect
cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed $dir0/Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8
cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed $dir0/Analysis/ATAC_seq/Union/R_output_Intersect/
cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed $dir0/Analysis/ATAC_seq/Union/R_output_FC75/
cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed $dir0/Analysis/ATAC_seq/Union/R_output_Close/NK_CD4_CD8/
cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed $dir0/Analysis/ATAC_seq/Union/R_output_Common/NK_CD4/
cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed $dir0/Analysis/ATAC_seq/Union/R_output_Common/NK_CD8/
cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed $dir0/Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8_FC_NK/
cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed $dir0/Analysis/ATAC_seq/Union/R_output_Common/Combine_All/


#Return peaks that are unique to each cell type specifically with a fold change of 1.5
#---------------------Returns background to be used for Homer Motif Analysis--------------
mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
for i in `ls *unique_all.bed`; do
  cell=${i%*_unique_all.bed}
  mkdir -p ${cell}_Motif
  mkdir -p ${cell}_Motif_Merge
  bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed -b $i | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > ${cell}_subtract.bed
  bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b $i | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > ${cell}_subtract_merge.bed
  #---------------------CREATE 6 COLUMN BED FILE--------------------------------------------
  ##This generates a 6 column bedgraph that includes unique peak names and adds strand info used for homer
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' $i > ${cell}_unique_all_5_column.bed ##Adds a unique peak name and strand info
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' ${cell}_subtract.bed > ${cell}_subtract_5_column.bed
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${cell}_subtract_merge.bed > ${cell}_subtract_5_column_merge.bed
  #findMotifsGenome.pl ${cell}_unique_all_5_column.bed mm10 ${cell}_Motif -size given -p 35 -bg ${cell}_subtract_5_column.bed	##Does motif analysis at unique ATAC seq peaks
  #findMotifsGenome.pl ${cell}_unique_all_5_column.bed mm10 ${cell}_Motif_Merge -size given -p 35 -bg ${cell}_subtract_5_column_merge.bed	##Does motif analysis at unique ATAC seq peaks
done

#Return peaks that are unique to each cell type specifically with a fold change of 1
cd $dir0/Analysis/ATAC_seq/Union/R_output_FC1
mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
for i in `ls *unique_all_FC1.bed`; do
  cell=${i%*_unique_all_FC1.bed}
  mkdir -p ${cell}_Motif_FC1
  mkdir -p ${cell}_Motif_Merge_FC1
  bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed -b $i | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > ${cell}_FC1_subtract.bed
  bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b $i | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > ${cell}_FC1_subtract_merge.bed
  #---------------------CREATE 6 COLUMN BED FILE--------------------------------------------
  ##This generates a 6 column bedgraph that includes unique peak names and adds strand info used for homer
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' $i > ${cell}_unique_all_FC1_5_column.bed ##Adds a unique peak name and strand info
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' ${cell}_FC1_subtract.bed > ${cell}_FC1_subtract_5_column.bed
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${cell}_FC1_subtract_merge.bed > ${cell}_FC1_subtract_5_column_merge.bed
  #findMotifsGenome.pl ${cell}_unique_all_FC1_5_column.bed mm10 ${cell}_Motif_FC1 -size given -p 35 -bg ${cell}_FC1_subtract_5_column.bed	##Does motif analysis at unique ATAC seq peaks
  #findMotifsGenome.pl ${cell}_unique_all_FC1_5_column.bed mm10 ${cell}_Motif_Merge_FC1 -size given -p 35 -bg ${cell}_FC1_subtract_5_column_merge.bed	##Does motif analysis at unique ATAC seq peaks

done

#Return peaks that are unique to each cell type specifically with a fold change of 0.75
cd $dir0/Analysis/ATAC_seq/Union/R_output_FC75
mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
for i in `ls *unique_all_FC75.bed`; do
  cell=${i%*_unique_all_FC75.bed}
  mkdir -p ${cell}_Motif_Merge_FC75
  bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b $i | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > ${cell}_FC75_subtract_merge.bed
  #---------------------CREATE 6 COLUMN BED FILE--------------------------------------------
  ##This generates a 6 column bedgraph that includes unique peak names and adds strand info used for homer
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' $i > ${cell}_unique_all_FC75_5_column.bed ##Adds a unique peak name and strand info
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${cell}_FC75_subtract_merge.bed > ${cell}_FC75_subtract_5_column_merge.bed
  #findMotifsGenome.pl ${cell}_unique_all_FC75_5_column.bed mm10 ${cell}_Motif_Merge_FC75 -size given -p 35 -bg ${cell}_FC75_subtract_5_column_merge.bed	##Does motif analysis at unique ATAC seq peaks
done

total_window_bp=2000
bin_size=20


#Return peaks that are common to CD4 and CD8, but unique to all others, but not NK cells
cd $dir0/Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8
mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b CD4_CD8.bed -wa > CD4_CD8_intersect.bed
mkdir -p CD4_CD8_Motif_Common
mkdir -p CD4_CD8_Motif_Common_Bit
bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b CD4_CD8_intersect.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > CD4_CD8_subtract_merge.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' CD4_CD8_subtract_merge.bed > CD4_CD8_subtract_merge_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' CD4_CD8_intersect.bed > CD4_CD8_intersect_5_column_merge.bed
#findMotifsGenome.pl CD4_CD8_intersect_5_column_merge.bed mm10 CD4_CD8_Motif_Common/ -size given -p 35 -bg CD4_CD8_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks
#findMotifsGenome.pl CD4_CD8_intersect_5_column_merge.bed mm10 CD4_CD8_Motif_Common_Bit/ -bits -size given -p 35 -bg CD4_CD8_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks


#annotatePeaks.pl CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m known3_ETS.motif known14_Runx1Jurkat.motif known_TCF3.motif > CD4_CD8_intersect_motif.count
#annotatePeaks.pl CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif motif1_TCF.motif > CD4_CD8_intersect_motif_new.count
#annotatePeaks.pl CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif Tcf1_ChIP_seq.motif > CD4_CD8_intersect_motif_new_chip.count

#Return peaks that are common to CD4 and CD8, but unique to all others, including NK cells
cd $dir0/Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8_FC_NK
mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b CD4_CD8.bed -wa > CD4_CD8_intersect.bed
mkdir -p CD4_CD8_FC_NK_Motif_Common
mkdir -p CD4_CD8_FC_NK_Motif_Common_Bit
bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b CD4_CD8_intersect.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > CD4_CD8_subtract_merge.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' CD4_CD8_subtract_merge.bed > CD4_CD8_subtract_merge_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' CD4_CD8_intersect.bed > CD4_CD8_intersect_5_column_merge.bed
#findMotifsGenome.pl CD4_CD8_intersect_5_column_merge.bed mm10 CD4_CD8_FC_NK_Motif_Common -size given -p 35 -bg CD4_CD8_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks
#findMotifsGenome.pl CD4_CD8_intersect_5_column_merge.bed mm10 CD4_CD8_FC_NK_Motif_Common_Bit -bits -size given -p 35 -bg CD4_CD8_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks

#annotatePeaks.pl CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m known3_ETS.motif known14_Runx1Jurkat.motif known_TCF3.motif > CD4_CD8_intersect_motif.count
#annotatePeaks.pl CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif motif1_TCF.motif > CD4_CD8_intersect_motif_new.count
#annotatePeaks.pl CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif Tcf1_ChIP_seq.motif > CD4_CD8_intersect_motif_new_chip.count


#Return peaks that are common to NK and CD4, but unique to all others
cd $dir0/Analysis/ATAC_seq/Union/R_output_Common/NK_CD4
mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b NK_CD4.bed -wa > NK_CD4_intersect.bed
mkdir -p NK_CD4_Motif_Common
mkdir -p NK_CD4_Motif_Common_Bit
bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b NK_CD4_intersect.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > NK_CD4_subtract_merge.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' NK_CD4_subtract_merge.bed > NK_CD4_subtract_merge_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' NK_CD4_intersect.bed > NK_CD4_intersect_5_column_merge.bed
#findMotifsGenome.pl NK_CD4_intersect_5_column_merge.bed mm10 NK_CD4_Motif_Common/ -size given -p 35 -bg NK_CD4_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks
#findMotifsGenome.pl NK_CD4_intersect_5_column_merge.bed mm10 NK_CD4_Motif_Common_Bit/ -bits -size given -p 35 -bg NK_CD4_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks

#annotatePeaks.pl NK_CD4_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m known3_ETS.motif known14_Runx1Jurkat.motif known_TCF3.motif > NK_CD4_intersect_motif.count
#annotatePeaks.pl NK_CD4_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif motif1_TCF.motif > NK_CD4_intersect_motif_new.count
#annotatePeaks.pl NK_CD4_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif Tcf1_ChIP_seq.motif > NK_CD4_intersect_motif_new_chip.count

#Return peaks that are common to NK and CD8, but unique to all others
cd $dir0/Analysis/ATAC_seq/Union/R_output_Common/NK_CD8
mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b NK_CD8.bed -wa > NK_CD8_intersect.bed
mkdir -p NK_CD8_Motif_Common
mkdir -p NK_CD8_Motif_Common_Bit
bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b NK_CD8_intersect.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > NK_CD8_subtract_merge.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' NK_CD8_subtract_merge.bed > NK_CD8_subtract_merge_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' NK_CD8_intersect.bed > NK_CD8_intersect_5_column_merge.bed
#findMotifsGenome.pl NK_CD8_intersect_5_column_merge.bed mm10 NK_CD8_Motif_Common/ -size given -p 35 -bg NK_CD8_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks
#findMotifsGenome.pl NK_CD8_intersect_5_column_merge.bed mm10 NK_CD8_Motif_Common_Bit/ -bits -size given -p 35 -bg NK_CD8_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks

#annotatePeaks.pl NK_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m known3_ETS.motif known14_Runx1Jurkat.motif known_TCF3.motif > NK_CD8_intersect_motif.count
#annotatePeaks.pl NK_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif motif1_TCF.motif > NK_CD8_intersect_motif_new.count
#annotatePeaks.pl NK_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif Tcf1_ChIP_seq.motif > NK_CD8_intersect_motif_new_chip.count

#Return peaks that are common to NK, CD4, and CD8, but unique to all others
cd $dir0/Analysis/ATAC_seq/Union/R_output_Common/NK_CD4_CD8
mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b NK_CD4_CD8.bed -wa > NK_CD4_CD8_intersect.bed
mkdir -p NK_CD4_CD8_Motif_Common
mkdir -p NK_CD4_CD8_Motif_Common_Bit
bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b NK_CD4_CD8_intersect.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > NK_CD4_CD8_subtract_merge.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' NK_CD4_CD8_subtract_merge.bed > NK_CD4_CD8_subtract_merge_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' NK_CD4_CD8_intersect.bed > NK_CD4_CD8_intersect_5_column_merge.bed
#findMotifsGenome.pl NK_CD4_CD8_intersect_5_column_merge.bed mm10 NK_CD4_CD8_Motif_Common/ -size given -p 35 -bg NK_CD4_CD8_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks
#findMotifsGenome.pl NK_CD4_CD8_intersect_5_column_merge.bed mm10 NK_CD4_CD8_Motif_Common_Bit/ -bits -size given -p 35 -bg NK_CD4_CD8_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks

#annotatePeaks.pl NK_CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m known3_ETS.motif known14_Runx1Jurkat.motif known_TCF3.motif > NK_CD4_CD8_intersect_motif.count
#annotatePeaks.pl NK_CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif motif1_TCF.motif > NK_CD4_CD8_intersect_motif_new.count
#annotatePeaks.pl NK_CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif Tcf1_ChIP_seq.motif > NK_CD4_CD8_intersect_motif_new_chip.count

cd $dir0/Analysis/ATAC_seq/Union/R_output_Common
cat CD4_CD8/CD4_CD8.bed NK_CD4_CD8/NK_CD4_CD8.bed NK_CD4/NK_CD4.bed NK_CD8/NK_CD8.bed > Combine_All_CD4_CD8_NK_Peaks.bed
cp Combine_All_CD4_CD8_NK_Peaks.bed Combine_All/
cd $dir0/Analysis/ATAC_seq/Union/R_output_Common/Combine_All
mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b Combine_All_CD4_CD8_NK_Peaks.bed -wa > Combine_All_CD4_CD8_NK_Peaks_intersect.bed
mkdir -p Combine_Peaks_All_Motif
mkdir -p Combine_Peaks_All_Motif_Bit
bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b Combine_All_CD4_CD8_NK_Peaks_intersect.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > Combine_All_CD4_CD8_NK_Peaks_subtract_merge.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' Combine_All_CD4_CD8_NK_Peaks_subtract_merge.bed > Combine_All_CD4_CD8_NK_Peaks_subtract_merge_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' Combine_All_CD4_CD8_NK_Peaks_intersect.bed > Combine_All_CD4_CD8_NK_Peaks_intersect_5_column_merge.bed
#findMotifsGenome.pl Combine_All_CD4_CD8_NK_Peaks_intersect_5_column_merge.bed mm10 Combine_Peaks_All_Motif/ -size given -p 35 -bg Combine_All_CD4_CD8_NK_Peaks_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks
#findMotifsGenome.pl Combine_All_CD4_CD8_NK_Peaks_intersect_5_column_merge.bed mm10 Combine_Peaks_All_Motif_Bit/ -bits -size given -p 35 -bg Combine_All_CD4_CD8_NK_Peaks_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks
#annotatePeaks.pl NK_CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m known3_ETS.motif known14_Runx1Jurkat.motif known_TCF3.motif > NK_CD4_CD8_intersect_motif.count
#annotatePeaks.pl NK_CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif motif1_TCF.motif > NK_CD4_CD8_intersect_motif_new.count
#annotatePeaks.pl NK_CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif Tcf1_ChIP_seq.motif > NK_CD4_CD8_intersect_motif_new_chip.count


#A test one where the intersection of the unique peaks were taken to the merged, union peak file to return the entire peak instead of the fragment
cd $dir0/Analysis/ATAC_seq/Union/R_output_Intersect
mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
for i in `ls *unique_all_FC1.bed`; do
  cell=${i%*_unique_all_FC1.bed}
  bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b $i -wa > ${cell}_unique_merge_intersect.bed
  mkdir -p ${cell}_Motif_Intersect
  bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b ${cell}_unique_merge_intersect.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > ${cell}_subtract_merge.bed
  #---------------------CREATE 6 COLUMN BED FILE--------------------------------------------
  ##This generates a 6 column bedgraph that includes unique peak names and adds strand info used for homer
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${cell}_subtract_merge.bed > ${cell}_subtract_merge_5_column.bed ##Adds a unique peak name and strand info
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${cell}_unique_merge_intersect.bed > ${cell}_unique_merge_intersect_5_column_merge.bed
  #findMotifsGenome.pl ${cell}_unique_merge_intersect_5_column_merge.bed mm10 ${cell}_Motif_Intersect/ -size given -p 35 -bg ${cell}_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks
done

#A test one where peaks that were closed in NK, CD4, and CD8 but open in B, Lsk, and Mono
#Did not work
cd $dir0/Analysis/ATAC_seq/Union/R_output_Close/NK_CD4_CD8
total_window_bp=2000
bin_size=20

mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b Lsk_B_Mono_close.bed -wa > Lsk_B_Mono_close_intersect.bed
mkdir -p NK_CD4_CD8_Motif_Close
bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b Lsk_B_Mono_close_intersect.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > Lsk_B_Mono_close_intersect_subtract.bed
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' Lsk_B_Mono_close_intersect_subtract.bed > Lsk_B_Mono_close_intersect_subtract_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' Lsk_B_Mono_close_intersect.bed > Lsk_B_Mono_close_intersect_5_column_merge.bed
#findMotifsGenome.pl Lsk_B_Mono_close_intersect_5_column_merge.bed mm10 NK_CD4_CD8_Motif_Close/ -size given -p 35 -bg Lsk_B_Mono_close_intersect_subtract_5_column.bed	##Does motif analysis at unique ATAC seq peaks
#annotatePeaks.pl CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m known3_ETS.motif known14_Runx1Jurkat.motif known_TCF3.motif > CD4_CD8_intersect_motif.count
#annotatePeaks.pl CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif motif1_TCF.motif > CD4_CD8_intersect_motif_new.count
#annotatePeaks.pl CD4_CD8_intersect.bed mm10 -multi -size ${total_window_bp} -hist $bin_size -m ETS1_ETS_Jurkat-ETS1-ChIP-Seq_GSE17954.motif RUNX1_Runt_Jurkat-RUNX1-ChIP-Seq_GSE29180.motif Tcf1_ChIP_seq.motif > CD4_CD8_intersect_motif_new_chip.count

