#!/bin/bash
cd /mnt/data1/John/Pioneer_Factors/Analysis/ATAC_seq/Union/R_output_Common/Intersect
mkdir -p motif_search
cp /mnt/data1/John/Pioneer_Factors/Analysis/ATAC_seq/Union/R_output_Common/NK_CD4_CD8/Tcf1_ChIP_seq.motif .
#cp /mnt/data1/John/Pioneer_Factors/Analysis/ATAC_seq/Union/R_output_Common/NK_CD4_CD8/motif1_TCF.motif .
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' Common_Peaks_All.bed > Common_Peaks_All_5_column.bed ##Adds a unique peak name and strand info
findMotifsGenome.pl Common_Peaks_All_5_column.bed mm10 motif_search/ -size given -find Tcf1_ChIP_seq.motif > outputfile.txt
#findMotifsGenome.pl Common_Peaks_All_5_column.bed mm10 motif_search/ -size given -find motif1_TCF.motif > outputfile.txt
#annotatePeaks.pl Common_Peaks_All_5_column.bed mm10 -size given -m Tcf1_ChIP_seq.motif > outputfile_annotate.txt

awk 'BEGIN {OFS="\t"} {print $1}' outputfile.txt > Peaks_with_motif.txt ##Adds a unique peak name and strand info
grep -Fwf Peaks_with_motif.txt Common_Peaks_All_5_column.bed > Peaks_with_motif_intersect.bed
awk '!seen[$4] {print} {++seen[$4]}' Peaks_with_motif_intersect.bed > Peaks_with_motif_intersect_dedup.bed
cp $datadir/ChIP_seq/macs/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph .
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph > Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed ##Adds a unique peak name and strand info

bedtools intersect -a Peaks_with_motif_intersect.bed -b Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed > Common_Peaks_All_Tcf1.bed

