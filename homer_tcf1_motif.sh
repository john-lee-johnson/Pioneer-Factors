#!/bin/bash

#----------------TCF-1 MOTIF CHIP-SEQ-----------------------------------------------------
#Generate a TCF-1 consensus motif from the Tcf-1 ChIP-seq using Homer
cd $dir0/Analysis/ChIP_seq/Transcription_Factor/Motif/Tcf1

#Copy over the Tcf-1 ChIP-seq peak file
cp $datadir/ChIP_seq/macs/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph .
#Generate a 5 column Homer bed file
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph > Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed ##Adds a unique peak name and strand info
#Homer motif search
mkdir Homer_Motif_Tcf1_Bit
#findMotifsGenome.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 Homer_Motif_Tcf1 -size 200
findMotifsGenome.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 Homer_Motif_Tcf1_Bit -bits -size 200
#Copy the de novo motif to the folders used for ATAC-seq peak motif analysis
cp $dir0/Analysis/ChIP_seq/Transcription_Factor/Motif/Tcf1/Homer_Motif_Tcf1/homerResults/motif2.motif $dir0/Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8/Tcf1_ChIP_seq.motif
cp $dir0/Analysis/ChIP_seq/Transcription_Factor/Motif/Tcf1/Homer_Motif_Tcf1/homerResults/motif2.motif $dir0/Analysis/ATAC_seq/Union/R_output_Common/NK_CD4_CD8/Tcf1_ChIP_seq.motif

#-------------OVERLAP OF ATAC SEQ PEAKS AND TCF-1 BINDING---------------------------------
cd $dir0/Analysis/ATAC_seq/Union/R_output_Common/Intersect
#Copy over the Tcf-1 binding sites
cp $dir0/Analysis/ChIP_seq/Transcription_Factor/Motif/Tcf1/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed .
#Copy over the CD4 and CD8 common Peaks
cp ../CD4_CD8/CD4_CD8_intersect_5_column_merge.bed .
#Copy over the NK, CD4 and CD8 common Peaks
cp ../NK_CD4_CD8/NK_CD4_CD8_intersect_5_column_merge.bed .
#Merge the CD4 and CD8 common peaks with the NK, CD4, and CD8 common peaks
cat *5_column_merge.bed > Common_Peaks_All.bed
#Find the intersect of the merged common peak file with the Tcf1 ChIP-seq peaks
bedtools intersect -a Common_Peaks_All.bed -b Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed > Common_Peaks_All_Tcf1.bed
#Find the intersect of the CD4 and CD8 common peak file with the Tcf1 ChIP-seq peaks
#bedtools intersect -a CD4_CD8_intersect_5_column_merge.bed -b Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed > CD4_CD8_Tcf1.bed
#Find the intersect of the CD4, CD8 and NK common peak file with the Tcf1 ChIP-seq peaks
#bedtools intersect -a NK_CD4_CD8_intersect_5_column_merge.bed -b Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed > NK_CD4_CD8_Tcf1.bed
#Find the intersect of the merged common peak file (with a window of 2kb) with the Tcf1 ChIP-seq peaks
#bedtools window -w 1000 -a Common_Peaks_All.bed -b Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed > Common_Peaks_All_Tcf1_Window.bed
#Find the intersect of the CD4 and CD8 common peak file (with a window of 2kb) with the Tcf1 ChIP-seq peaks
#bedtools window -w 1000 -a CD4_CD8_intersect_5_column_merge.bed -b Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed > CD4_CD8_Tcf1_Window.bed
#Find the intersect of the CD4, CD8 and NK common peak file (with a window of 2kb) with the Tcf1 ChIP-seq peaks
#bedtools window -w 1000 -a NK_CD4_CD8_intersect_5_column_merge.bed -b Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed > NK_CD4_CD8_Tcf1_Window.bed

#---------------REDO MOTIF ANALYSIS WITH OVERLAP------------------------------------------
#With the intersect between the common peaks and the Tcf-1 ChIP-seq, perform motif analysis
#cd $dir0/Analysis/ATAC_seq/Union/R_output_Common/Intersect/Motif_Overlap
cp $dir0/Analysis/ATAC_seq/Union/R_output_Common/NK_CD4_CD8/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed .
bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b Common_Peaks_All_Tcf1.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > Common_Peaks_All_Tcf1_subtract_merge.bed

awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' Common_Peaks_All_Tcf1.bed > Common_Peaks_All_Tcf1_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' Common_Peaks_All_Tcf1_subtract_merge.bed > Common_Peaks_All_Tcf1_subtract_merge_5_column.bed
#findMotifsGenome.pl Common_Peaks_All_Tcf1_5_column.bed mm10 Motif_Overlap/ -size given -p 35 -bg Common_Peaks_All_Tcf1_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks
