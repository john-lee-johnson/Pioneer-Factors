#!/bin/bash
cd /mnt/data1/VahediLab/PTF_Team/Data/ChIP_seq/other
awk '{FS=",";OFS="\t"; print $1,$4,$5,$6}' GSE56393_ETS1_DN_mm9.csv > GSE56393_ETS1_DN_mm9_temp.bedgraph
tail -n +2 GSE56393_ETS1_DN_mm9_temp.bedgraph > GSE56393_ETS1_DN_mm9.bedgraph
awk '{FS=",";OFS="\t"; print $1,$4,$5,$6}' GSE56393_ETS1_DP_mm9.csv > GSE56393_ETS1_DP_mm9_temp.bedgraph
tail -n +2 GSE56393_ETS1_DP_mm9_temp.bedgraph > GSE56393_ETS1_DP_mm9.bedgraph
rm *temp.bedgraph
CrossMap.py bed /mnt/data1/VahediLab/PTF_Team/Data/mm9ToMm10.over.chain.gz GSE56393_ETS1_DP_mm9.bedgraph GSE56393_ETS1_DP_mm10_unfixed.bedgraph
CrossMap.py bed /mnt/data1/VahediLab/PTF_Team/Data/mm9ToMm10.over.chain.gz GSE56393_ETS1_DN_mm9.bedgraph GSE56393_ETS1_DN_mm10_unfixed.bedgraph
awk 'OFS="\t"{if ($3-$2<0) {start=$3; end=$2} else {start=$2; end=$3} print $1,start,end,$4}' GSE56393_ETS1_DP_mm10_unfixed.bedgraph > DP_ChIP_seq_ETS1_Cauchy_mm10_GSE56393.bedgraph
#awk 'OFS="\t"{if ($3-$2==0) {start=$2-1; end=$3} else {start=$2; end=$3} print $1,start,end,$4}' GSE56393_ETS1_DP_mm10_temp.bedgraph > DP_ChIP_seq_ETS1_Cauchy_mm10_GSE56393.bedgraph

awk 'OFS="\t"{if ($3-$2<0) {start=$3; end=$2} else {start=$2; end=$3} print $1,start,end,$4}' GSE56393_ETS1_DN_mm10_unfixed.bedgraph > DN_ChIP_seq_ETS1_Cauchy_mm10_GSE56393.bedgraph
#awk 'OFS="\t"{if ($3-$2==0) {start=$2-1; end=$3} else {start=$2; end=$3} print $1,start,end,$4}' GSE56393_ETS1_DN_mm10_temp.bedgraph > DN_ChIP_seq_ETS1_Cauchy_mm10_GSE56393.bedgraph

bedtools unionbedg -i DN_ChIP_seq_ETS1_Cauchy_mm10_GSE56393.bedgraph DP_ChIP_seq_ETS1_Cauchy_mm10_GSE56393.bedgraph > ETS1_combined_DN_DP_mm10.bedgraph
awk '{FS="\t";OFS="\t"; sum=$4+$5; print $1,$2,$3,sum}' ETS1_combined_DN_DP_mm10.bedgraph > ETS1_combined_DN_DP_mm10_sum_temp.bedgraph
awk 'OFS="\t"{if ($3-$2<0) {start=$3; end=$2} else {start=$2; end=$3} print $1,start,end,$4}' ETS1_combined_DN_DP_mm10_sum_temp.bedgraph > ETS1_combined_DN_DP_mm10_sum.bedgraph

#awk 'OFS="\t"{if ($3-$2<0) {start=$3; end=$2} else {start=$2; end=$3} print $1,start,end,$4}' ETS1_combined_DN_DP_mm10_sum_temp.bedgraph > ETS1_combined_DN_DP_mm10_sum.bedgraph
sort -k1,1 -k2,2n ETS1_combined_DN_DP_mm10_sum.bedgraph > ETS1_combined_DN_DP_mm10_sorted.bedgraph
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' ETS1_combined_DN_DP_mm10_sorted.bedgraph > ETS1_combined_DN_DP_mm10_sorted_5_column.bed ##Adds a unique peak name and strand info
bedtools merge -i ETS1_combined_DN_DP_mm10_sorted_5_column.bed > ETS1_combined_DN_DP_mm10_sorted_merged_5_column.bed

rm GSE56393_ETS1_DN_mm10_unfixed.bedgraph
rm GSE56393_ETS1_DP_mm10_unfixed.bedgraph
cp DP_ChIP_seq_ETS1_Cauchy_mm10_GSE56393.bedgraph DP_ChIP_seq_ETS1_Cauchy_mm10_GSE56393_tag_peaks.bedgraph
cp DN_ChIP_seq_ETS1_Cauchy_mm10_GSE56393.bedgraph DN_ChIP_seq_ETS1_Cauchy_mm10_GSE56393_tag_peaks.bedgraph
#chr10	3516739	3511739	peak744	0	+

dir1=/mnt/data1/John/Pioneer_Factors/Analysis/ATAC_seq/Union/R_output_Common_1.5
dir2=/mnt/data1/VahediLab/PTF_Team/Data/ChIP_seq

rm -f $dir1/tcf1_ets1.txt
cd $dir1
for i in `ls -d */`; do
#echo "$i"
filename=$(echo $i | cut -d '/' -f1)
echo "$filename $(cat $dir1/$filename/${filename}.bed | wc -l) ETS1_DN $(bedtools window -a $dir1/$filename/${filename}.bed -b $dir2/other/DN_ChIP_seq_ETS1_Cauchy_mm10_GSE56393.bedgraph -w 2 | wc -l) ETS1_DP $(bedtools window -a $dir1/$filename/${filename}.bed -b $dir2/other/DP_ChIP_seq_ETS1_Cauchy_mm10_GSE56393.bedgraph -w 2 | wc -l) TCF1 $(bedtools window -a $dir1/$filename/${filename}.bed -b $dir2/macs/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph -w 2 | wc -l)" >> $dir1/tcf1_ets1.txt

done

dir1=/mnt/data1/John/Pioneer_Factors/Analysis/ATAC_seq/Union/R_output_Common_Merge1.5
rm -f $dir1/tcf1_ets1.txt
cd $dir1
for i in `ls -d */`; do
#echo "$i"
filename=$(echo $i | cut -d '/' -f1)
echo "$filename $(cat $dir1/$filename/${filename}.bed | wc -l) ETS1_DN $(bedtools window -a $dir1/$filename/${filename}.bed -b $dir2/other/DN_ChIP_seq_ETS1_Cauchy_mm10_GSE56393.bedgraph -w 2 | wc -l) ETS1_DP $(bedtools window -a $dir1/$filename/${filename}.bed -b $dir2/other/DP_ChIP_seq_ETS1_Cauchy_mm10_GSE56393.bedgraph -w 2 | wc -l) TCF1 $(bedtools window -a $dir1/$filename/${filename}.bed -b $dir2/macs/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph -w 2 | wc -l)" >> $dir1/tcf1_ets1.txt

done
#/ChIP_seq/macs/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph

