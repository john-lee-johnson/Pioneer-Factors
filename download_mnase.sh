#!/bin/bash


#Set Directories
infodir=/mnt/data1/John/Pioneer-Factors/info
maindir=/mnt/data1/John/Pioneer_Factors
scriptdir=/mnt/data1/John/Pioneer-Factors
logdir=/mnt/data1/John/Pioneer-Factors/logs
filedir=/mnt/data1/John/Pioneer-Factors/sample_files
paralleldir=/mnt/data1/John/Pioneer-Factors/parallel_commands
r_dir=/mnt/data1/John/Pioneer-Factors/R_Scripts
datadir=/mnt/data1/VahediLab/PTF_Team/Data
datadir2=/mnt/data1/VahediLab/PTF_Team/Data2
dir0=$maindir

#Set datasets
tcf1_eml="EML_ChIP_seq_Tcf1_Wu_mm10_GSM773994"
tcf1_thy="Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644"
ets1_dn="DN_ChIP_seq_ETS1_Cauchy_mm10_GSE56393"
ets1_dp="DP_ChIP_seq_ETS1_Cauchy_mm10_GSE56393"
ets1_union="ETS1_combined_DN_DP_mm10_sorted_merged_5_column"
runx1_dn="DN_ChIP_seq_Runx1_Cauchy_mm10_GSM1360735"
runx1_dp="DP_ChIP_seq_Runx1_Bergon_mm10_GSM1095815"
gata3_dn="DN_ChIP_seq_Gata3_Zhao_mm10_GSM523221"
gata3_dp_zhao="DP_ChIP_seq_Gata3_Zhao_mm10_GSM523222"
gata3_dp_roth="DP_ChIP_seq_Gata3_Rothenberg_mm10_GSM774297"

#Copy files
#Tcf1
cp $datadir2/ChIP_seq/macs/${tcf1_eml}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DP
cp $datadir/ChIP_seq/macs/${tcf1_thy}_macs_peaks_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DP
#Ets1
cp $datadir/ChIP_seq/other/${ets1_dn}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DP
cp $datadir/ChIP_seq/other/${ets1_dp}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DP
cp $datadir/ChIP_seq/other/${ets1_union}.bed $dir0/Analysis/MNase_seq/DP
#Runx1
cp $datadir2/ChIP_seq/macs/${runx1_dn}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DP
cp $datadir2/ChIP_seq/macs/${runx1_dp}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DP
#Gata3
cp $datadir2/ChIP_seq/macs/${gata3_dn}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DP
cp $datadir2/ChIP_seq/macs/${gata3_dp_zhao}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DP
cp $datadir2/ChIP_seq/macs//${gata3_dp_roth}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DP

#Copy files
#Tcf1
cp $datadir2/ChIP_seq/macs/${tcf1_eml}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DN
cp $datadir/ChIP_seq/macs/${tcf1_thy}_macs_peaks_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DN
#Ets1
cp $datadir/ChIP_seq/other/${ets1_dn}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DN
cp $datadir/ChIP_seq/other/${ets1_dp}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DN
cp $datadir/ChIP_seq/other/${ets1_union}.bed $dir0/Analysis/MNase_seq/DN
#Runx1
cp $datadir2/ChIP_seq/macs/${runx1_dn}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DN
cp $datadir2/ChIP_seq/macs/${runx1_dp}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DN
#Gata3
cp $datadir2/ChIP_seq/macs/${gata3_dn}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DN
cp $datadir2/ChIP_seq/macs/${gata3_dp_zhao}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DN
cp $datadir2/ChIP_seq/macs//${gata3_dp_roth}_tag_peaks.bedgraph $dir0/Analysis/MNase_seq/DN

cd $dir0/Analysis/Homer/Tag_Directories/MNase
#makeTagDirectory DN_MNase_seq_mm10_GSM1359851_shift "/mnt/data1/John/Pioneer_Factors/Data/MNase_seq/DN/Cauchy/DN_MNase_seq_Cauchy_mm10_GSM1359851_shift.bam" -genome mm10 -format sam  
#makeTagDirectory DP_MNase_seq_mm10_GSM1359852_shift "/mnt/data1/John/Pioneer_Factors/Data/MNase_seq/DP/Cauchy/DP_MNase_seq_Cauchy_mm10_GSM1359852_shift.bam" -genome mm10 -format sam
#makeTagDirectory DN_MNase_seq_mm10_GSM1359851_shift "/mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851.bam" -genome mm10 -format sam  
#makeTagDirectory DP_MNase_seq_mm10_GSM1359852_shift "/mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852.bam" -genome mm10 -format sam

cd $dir0/Analysis/MNase_seq/DP/
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' ${tcf1_eml}_tag_peaks.bedgraph > ${tcf1_eml}_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' ${tcf1_thy}_macs_peaks_tag_peaks.bedgraph > ${tcf1_thy}_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' ${ets1_dn}_tag_peaks.bedgraph > ${ets1_dn}_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' ${ets1_dp}_tag_peaks.bedgraph > ${ets1_dp}_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' ${runx1_dn}_tag_peaks.bedgraph > ${runx1_dn}_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' ${runx1_dp}_tag_peaks.bedgraph > ${runx1_dp}_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' ${gata3_dn}_tag_peaks.bedgraph > ${gata3_dn}_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' ${gata3_dp_zhao}_tag_peaks.bedgraph > ${gata3_dp_zhao}_5_column.bed ##Adds a unique peak name and strand info
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' ${gata3_dp_roth}_tag_peaks.bedgraph > ${gata3_dp_roth}_5_column.bed ##Adds a unique peak name and strand info


cd $dir0/Analysis/MNase_seq/DP
#python /mnt/data1/bin/danpos-2.2.2/danpos.py dpos /mnt/data1/John/Pioneer_Factors/Data/MNase_seq/DP/Cauchy/DP_MNase_seq_Cauchy_mm10_GSM1359852.bam


#python /mnt/data1/bin/danpos-2.2.2/danpos.py profile mnt_data1_John_Pioneer_Factors_Data_MNase_seq_DP_Cauchy_DP_MNase_seq_Cauchy_mm10_GSM1359852.smooth.wig --bed3file_paths Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed --name Tcf1_Thy

cd $dir0/Analysis/MNase_seq/DN
python /mnt/data1/bin/danpos-2.2.2/danpos.py dpos /mnt/data1/John/Pioneer_Factors/Data/MNase_seq/DN/Cauchy/DN_MNase_seq_Cauchy_mm10_GSM1359851.bam

#-----------Normalizing to total number of reads----------------------
#annotatePeaks.pl ${tcf1_thy}_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Tcf1_Thy_histfile_reads.txt
#annotatePeaks.pl ${tcf1_eml}_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Tcf1_EML_histfile_reads.txt
#annotatePeaks.pl ${ets1_union}.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Ets1_Union_histfile_reads.txt
#annotatePeaks.pl ${ets1_dn}_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Ets1_DN_histfile_reads.txt
#annotatePeaks.pl ${ets1_dp}_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Ets1_DP_histfile_reads.txt
#annotatePeaks.pl ${runx1_dn}_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Runx1_DN_histfile_reads.txt
#annotatePeaks.pl ${runx1_dp}_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Runx1_DP_histfile_reads.txt

#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Ets1_histfile_reads.txt
#annotatePeaks.pl DP_ChIP_seq_Runx1_Bergon_mm10_GSM1095815_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Runx1_histfile_reads.txt
#annotatePeaks.pl DP_ChIP_seq_Gata3_Rothenberg_mm10_GSM774297_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Gata3_histfile_reads.txt
#annotatePeaks.pl EML_ChIP_seq_Tcf1_Wu_mm10_GSM773994_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > EML_histfile_reads.txt

#-----------Normalizing to total number of reads shifted bam----------
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > Tcf1_histfile_shift_reads.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > Ets1_histfile_shift_reads.txt
#annotatePeaks.pl DP_ChIP_seq_Runx1_Bergon_mm10_GSM1095815_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > Runx1_histfile_shift_reads.txt
#annotatePeaks.pl DP_ChIP_seq_Gata3_Rothenberg_mm10_GSM774297_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > Gata3_histfile_shift_reads.txt
#annotatePeaks.pl EML_ChIP_seq_Tcf1_Wu_mm10_GSM773994_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > EML_histfile_shift_reads.txt

#-----------Normalizing to 10E7----------------------
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -norm 1e7 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852 > Tcf1_histfile_norm.txt

#-----------Normalizing to local area raw bam----------------
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Tcf1_histfile_hist.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Ets1_histfile_hist.txt
#annotatePeaks.pl DP_ChIP_seq_Runx1_Bergon_mm10_GSM1095815_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Runx1_histfile_hist.txt
#annotatePeaks.pl DP_ChIP_seq_Gata3_Rothenberg_mm10_GSM774297_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Gata3_histfile_hist.txt
#annotatePeaks.pl EML_ChIP_seq_Tcf1_Wu_mm10_GSM773994_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > EML_histfile_hist.txt

#-----------Normalizing to local area shifted bam----------------
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > Tcf1_histfile_shift_hist.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > Ets1_histfile_shift_hist.txt
#annotatePeaks.pl DP_ChIP_seq_Runx1_Bergon_mm10_GSM1095815_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > Runx1_histfile_shift_hist.txt
#annotatePeaks.pl DP_ChIP_seq_Gata3_Rothenberg_mm10_GSM774297_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > Gata3_histfile_shift_hist.txt
#annotatePeaks.pl EML_ChIP_seq_Tcf1_Wu_mm10_GSM773994_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > EML_histfile_shift_hist.txt


cd $dir0/Analysis/MNase_seq/DP
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > Ets1_histfile_reads.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -norm 60981718 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DP_MNase_seq_Cauchy_mm10_GSM1359852 > Ets1_histfile_reads.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -norm 1e7 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852 > Ets1_histfile_norm.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -histNorm 10 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852 > Ets1_histfile_hist.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -norm 68409187 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > Ets1_histfile_shift.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -histNorm 10 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DP_MNase_seq_mm10_GSM1359852_shift > Ets1_histfile_shift.txt



cd $dir0/Analysis/MNase_seq/DN

#Tcf1
#cp $tcf1_thy .
#cp $tcf1_eml .
#Ets1
#cp $ets1_dn .
#cp $ets1_dp .
#cp $ets1_union .
#Runx1
#cp $runx1_dn .
#cp $runx1_dn .
#Gata3
#cp $gata3_dn .
#cp $gata3_dp_zhao .
#cp $gata3_dp_roth .

#Tcf1
echo "$tcf1_thy"
echo "$tcf1_eml"
#Ets1
echo "$ets1_dn"
echo "$ets1_dp"
echo "$ets1_union"
#Runx1
echo "$runx1_dn"
echo "$runx1_dn"
#Gata3
echo "$gata3_dn"
echo "$gata3_dp_zhao"
echo "$gata3_dp_roth"

#awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph > Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed ##Adds a unique peak name and strand info
#awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' DP_ChIP_seq_Runx1_Bergon_mm10_GSM1095815_tag_peaks.bedgraph > DP_ChIP_seq_Runx1_Bergon_mm10_GSM1095815_5_column.bed ##Adds a unique peak name and strand info
#awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' DP_ChIP_seq_Gata3_Rothenberg_mm10_GSM774297_tag_peaks.bedgraph > DP_ChIP_seq_Gata3_Rothenberg_mm10_GSM774297_5_column.bed ##Adds a unique peak name and strand info
#awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' EML_ChIP_seq_Tcf1_Wu_mm10_GSM773994_tag_peaks.bedgraph > EML_ChIP_seq_Tcf1_Wu_mm10_GSM773994_5_column.bed ##Adds a unique peak name and strand info

#-----------Normalizing to total number of reads----------------------
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -norm 68409187  -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851 > Tcf1_histfile_reads.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -norm 68409187  -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851 > Ets1_histfile_reads.txt
#annotatePeaks.pl DP_ChIP_seq_Runx1_Bergon_mm10_GSM1095815_5_column.bed mm10 -norm 68409187  -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851 > Runx1_histfile_reads.txt
#annotatePeaks.pl DP_ChIP_seq_Gata3_Rothenberg_mm10_GSM774297_5_column.bed mm10 -norm 68409187  -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851 > Gata3_histfile_reads.txt
#annotatePeaks.pl EML_ChIP_seq_Tcf1_Wu_mm10_GSM773994_5_column.bed mm10 -norm 68409187  -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851 > EML_histfile_reads.txt

#-----------Normalizing to total number of reads shifted bam----------
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -norm 68409187 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Tcf1_histfile_shift_reads.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -norm 68409187 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Ets1_histfile_shift_reads.txt
#annotatePeaks.pl DP_ChIP_seq_Runx1_Bergon_mm10_GSM1095815_5_column.bed mm10 -norm 68409187 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Runx1_histfile_shift_reads.txt
#annotatePeaks.pl DP_ChIP_seq_Gata3_Rothenberg_mm10_GSM774297_5_column.bed mm10 -norm 68409187 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Gata3_histfile_shift_reads.txt
#annotatePeaks.pl EML_ChIP_seq_Tcf1_Wu_mm10_GSM773994_5_column.bed mm10 -norm 68409187 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > EML_histfile_shift_reads.txt


#-----------Normalizing to local area raw bam----------------
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851 > Tcf1_histfile_hist.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851 > Ets1_histfile_hist.txt
#annotatePeaks.pl DP_ChIP_seq_Runx1_Bergon_mm10_GSM1095815_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851 > Runx1_histfile_hist.txt
#annotatePeaks.pl DP_ChIP_seq_Gata3_Rothenberg_mm10_GSM774297_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851 > Gata3_histfile_hist.txt
#annotatePeaks.pl EML_ChIP_seq_Tcf1_Wu_mm10_GSM773994_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851 > EML_histfile_hist.txt


#-----------Normalizing to local area shifted bam----------------
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Tcf1_histfile_shift_hist.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Ets1_histfile_shift_hist.txt
#annotatePeaks.pl DP_ChIP_seq_Runx1_Bergon_mm10_GSM1095815_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Runx1_histfile_shift_hist.txt
#annotatePeaks.pl DP_ChIP_seq_Gata3_Rothenberg_mm10_GSM774297_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Gata3_histfile_shift_hist.txt
#annotatePeaks.pl EML_ChIP_seq_Tcf1_Wu_mm10_GSM773994_5_column.bed mm10 -histNorm 5 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > EML_histfile_shift_hist.txt



#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -norm 68409187 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Tcf1_histfile_reads.txt
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -norm 68409187 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851 > Tcf1_histfile_reads.txt
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -norm 1e7 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851 > Tcf1_histfile_norm.txt
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -histNorm 10 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851 > Tcf1_histfile_hist.txt
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -norm 68409187 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Tcf1_histfile_shift.txt
#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -histNorm 10 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Tcf1_histfile_shift.txt


cd $dir0/Analysis/MNase_seq/DN
#cp $datadir/ChIP_seq/other/ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed .
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -norm 68409187 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Ets1_histfile_reads.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -norm 68409187 -size 4000 -hist 10 -d /mnt/data1/John/Pioneer_Factors/Analysis/Homer/Tag_Directories/GSE56360/DN_MNase_seq_Cauchy_mm10_GSM1359851 > Ets1_histfile_reads.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -norm 1e7 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851 > Ets1_histfile_norm.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -histNorm 10 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851 > Ets1_histfile_hist.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -norm 68409187 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Ets1_histfile_shift.txt
#annotatePeaks.pl ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed mm10 -histNorm 10 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DN_MNase_seq_mm10_GSM1359851_shift > Ets1_histfile_shift.txt
