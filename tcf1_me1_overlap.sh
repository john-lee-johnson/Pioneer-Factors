#!/bin/bash

#Intersect of H3K4me1 and Tcf1
mkdir -p /mnt/data1/John/pioneer/analysis/chip_seq/tcf1/intersect
cd /mnt/data1/John/pioneer/analysis/chip_seq/tcf1/intersect
dir0=/mnt/data1/John/pioneer/analysis/chip_seq/tcf1/macs
dir1=/mnt/data1/John/pioneer/analysis/chip_seq/macs
#tcf=/mnt/data1/John/pioneer/analysis/chip_seq/tcf1/macs/Thy_Tcf1_macs_peaks.bed
tcf=/mnt/data1/John/pioneer/analysis/chip_seq/tcf1/tag/Thy_Tcf1_window.bed
#tcf=/mnt/data1/John/pioneer/analysis/chip_seq/tcf1/tag/Thy_Tcf1_window_500.bed

bedtools intersect -a $tcf -b $dir1/LTHSC_H3K4me1_macs_14_peaks.bed > intersect_TCF1_LTHSC_H3K4me1.bed
bedtools intersect -a $tcf -b $dir1/LTHSC_H3K4me2_macs_14_peaks.bed > intersect_TCF1_LTHSC_H3K4me2.bed
bedtools intersect -a $tcf -b $dir1/STHSC_H3K4me1_macs_14_peaks.bed > intersect_TCF1_STHSC_H3K4me1.bed
bedtools intersect -a $tcf -b $dir1/STHSC_H3K4me2_macs_14_peaks.bed > intersect_TCF1_STHSC_H3K4me2.bed
bedtools intersect -a $tcf -b $dir1/MPP_H3K4me1_macs_14_peaks.bed > intersect_TCF1_MPP_H3K4me1.bed
bedtools intersect -a $tcf -b $dir1/MPP_H3K4me2_macs_14_peaks.bed > intersect_TCF1_MPP_H3K4me2.bed
bedtools intersect -a $tcf -b $dir1/CLP_H3K4me1_macs_14_peaks.bed > intersect_TCF1_CLP_H3K4me1.bed
bedtools intersect -a $tcf -b $dir1/CLP_H3K4me2_macs_14_peaks.bed > intersect_TCF1_CLP_H3K4me2.bed
bedtools intersect -a $tcf -b $dir1/CD4_H3K4me1_macs_14_peaks.bed > intersect_TCF1_CD4_H3K4me1.bed
bedtools intersect -a $tcf -b $dir1/CD4_H3K4me2_macs_14_peaks.bed > intersect_TCF1_CD4_H3K4me2.bed
bedtools intersect -a $tcf -b $dir1/CD8_H3K4me1_macs_14_peaks.bed > intersect_TCF1_CD8_H3K4me1.bed
bedtools intersect -a $tcf -b $dir1/CD8_H3K4me2_macs_14_peaks.bed > intersect_TCF1_CD8_H3K4me2.bed
bedtools intersect -a $tcf -b $dir1/LTHSC_H3K27Ac_macs_14_peaks.bed > intersect_TCF1_LTHSC_H3K27Ac.bed
bedtools intersect -a $tcf -b $dir1/STHSC_H3K27Ac_macs_14_peaks.bed > intersect_TCF1_STHSC_H3K27Ac.bed
bedtools intersect -a $tcf -b $dir1/MPP_H3K27Ac_macs_14_peaks.bed > intersect_TCF1_MPP_H3K27Ac.bed
bedtools intersect -a $tcf -b $dir1/CLP_H3K27Ac_macs_14_peaks.bed > intersect_TCF1_CLP_H3K27Ac.bed
bedtools intersect -a $tcf -b $dir1/CD4_H3K27Ac_macs_14_peaks.bed > intersect_TCF1_CD4_H3K27Ac.bed
bedtools intersect -a $tcf -b $dir1/CD8_H3K27Ac_macs_14_peaks.bed > intersect_TCF1_CD8_H3K27Ac.bed

