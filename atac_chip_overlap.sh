#!/bin/bash
#The purpose of this script is to look at the proportion of unique atac-seq peaks (over 1.5 fc) bound by a certain transcription factor
#The first step is to go into a new folder (Pioneer_Factors/Analysis/ATAC_seq/Overlap

#--------------Sets the working directories-----------------------------------------------
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

cd $dir0/Analysis/ATAC_seq/Overlap
rm -f $infodir/chip_seq_list.txt
#Next we set the ChIP-seq datasets we want to check
rm -f atac_chip_overlap.txt
rm -f atac_chip_overlap_clusters_0.txt
rm -f atac_chip_overlap_clusters_0_merge.txt
ets1_dp="/mnt/data1/VahediLab/PTF_Team/Data/ChIP_seq/other/DP_ChIP_seq_ETS1_Cauchy_mm10_GSE56393_tag_peaks.bedgraph"
ets1_dn="/mnt/data1/VahediLab/PTF_Team/Data/ChIP_seq/other/DN_ChIP_seq_ETS1_Cauchy_mm10_GSE56393_tag_peaks.bedgraph"
tcf1_dp="/mnt/data1/VahediLab/PTF_Team/Data/ChIP_seq/macs/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_tag_peaks.bedgraph"
tcf1_eml="/mnt/data1/VahediLab/PTF_Team/Data2/ChIP_seq/macs/EML_ChIP_seq_Tcf1_Wu_mm10_GSM773994_tag_peaks.bedgraph"
gata3_dp="/mnt/data1/VahediLab/PTF_Team/Data2/ChIP_seq/macs/DP_ChIP_seq_Gata3_Zhao_mm10_GSM523222_tag_peaks.bedgraph"
gata3_dn="/mnt/data1/VahediLab/PTF_Team/Data2/ChIP_seq/macs/DN_ChIP_seq_Gata3_Zhao_mm10_GSM523221_tag_peaks.bedgraph"
runx1_dp="/mnt/data1/VahediLab/PTF_Team/Data2/ChIP_seq/macs/DP_ChIP_seq_Runx1_Bergon_mm10_GSM1095815_tag_peaks.bedgraph"
runx1_dn="/mnt/data1/VahediLab/PTF_Team/Data2/ChIP_seq/macs/DN_ChIP_seq_Runx1_Cauchy_mm10_GSM1360735_tag_peaks.bedgraph"
pu1_dp="/mnt/data1/VahediLab/PTF_Team/Data2/ChIP_seq/macs/DP_ChIP_seq_PU.1_Rothenberg_mm10_GSM774294_tag_peaks.bedgraph"
pu1_dn1="/mnt/data1/VahediLab/PTF_Team/Data2/ChIP_seq/macs/FLDN1_ChIP_seq_PU.1_Rothenberg_mm10_GSM774291_tag_peaks.bedgraph"
pu1_dn2a="/mnt/data1/VahediLab/PTF_Team/Data2/ChIP_seq/macs/FLDN2a_ChIP_seq_PU.1_Rothenberg_mm10_GSM774292_tag_peaks.bedgraph"
pu1_dn2b="/mnt/data1/VahediLab/PTF_Team/Data2/ChIP_seq/macs/FLDN2b_ChIP_seq_PU.1_Rothenberg_mm10_GSM774293_tag_peaks.bedgraph"

echo "$ets1_dp" >> $infodir/chip_seq_list.txt
echo "$ets1_dn" >> $infodir/chip_seq_list.txt
echo "$gata3_dn" >> $infodir/chip_seq_list.txt
echo "$gata3_dp" >> $infodir/chip_seq_list.txt
echo "$runx1_dp" >> $infodir/chip_seq_list.txt
echo "$runx1_dn" >> $infodir/chip_seq_list.txt
echo "$pu1_dp" >> $infodir/chip_seq_list.txt
echo "$pu1_dn1" >> $infodir/chip_seq_list.txt
echo "$pu1_dn2a" >> $infodir/chip_seq_list.txt
echo "$pu1_dn2b" >> $infodir/chip_seq_list.txt
echo "$tcf1_eml" >> $infodir/chip_seq_list.txt
echo "$tcf1_dp" >> $infodir/chip_seq_list.txt

lines=$(wc -l $infodir/chip_seq_list.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
  line=$(sed -n "${i}p" < $infodir/chip_seq_list.txt)
  peakFile=$(echo "$line")
  fullName=$(echo "$peakFile" | rev | cut -d'/' -f1 | rev)
  filename=${fullName%_tag_peaks.bedgraph}
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' $peakFile > ${filename}_5_column.bed ##Adds a unique peak name and strand info
done

#for a in `ls -d $dir0/Analysis/ATAC_seq/Union/R_output_Common_1.5/*/ | grep -v Clusters`; do
#  location=$(echo "$a" | rev | cut -d'/' -f2 | rev)
#  lines=$(wc -l $infodir/chip_seq_list.txt | cut -d' ' -f1)
#  for ((i=1; i<=$lines; i++)); do
#    line=$(sed -n "${i}p" < $infodir/chip_seq_list.txt)
#    peakFile=$(echo "$line")
#    fullName=$(echo "$peakFile" | rev | cut -d'/' -f1 | rev)
#    filename=${fullName%_tag_peaks.bedgraph}
#    cell=$(echo "$filename" | cut -d'_' -f1)
#    chip=$(echo "$filename" | cut -d'_' -f4)
#    investigator=$(echo "$filename" | cut -d'_' -f5)
#    bound=$(bedtools window -a $a/${location}.bed -b ${filename}_5_column.bed -w 2 | wc -l)
#    total=$(cat $a/${location}.bed | wc -l)
#    results=$( echo "$bound / $total" | bc -l )
#    echo "$location $cell $chip $investigator $bound $total $results"
#    echo "$location $cell $chip $investigator $bound $total $results" >> atac_chip_overlap.txt
#  done
#done

#for a in `ls $dir0/Analysis/ATAC_seq/Union/R_output_Common_0/Clusters/*15-clusters.bed`; do
#  location=$(echo "$a" | rev | cut -d'/' -f2 | rev)
#  cluster=$(echo "$a" | rev | cut -d'/' -f1 | rev | cut -d'.' -f3)
#  clusterTotal=$(echo "$a" | rev | cut -d'/' -f1 | rev | cut -d'.' -f5 | cut -d'-' -f1)
#  lines=$(wc -l $infodir/chip_seq_list.txt | cut -d' ' -f1)
#  for ((i=1; i<=$lines; i++)); do
#    line=$(sed -n "${i}p" < $infodir/chip_seq_list.txt)
#    peakFile=$(echo "$line")
#    fullName=$(echo "$peakFile" | rev | cut -d'/' -f1 | rev)
#    filename=${fullName%_tag_peaks.bedgraph}
#    cell=$(echo "$filename" | cut -d'_' -f1)
#    chip=$(echo "$filename" | cut -d'_' -f4)
#    investigator=$(echo "$filename" | cut -d'_' -f5)
#    bound=$(bedtools window -a $a -b ${filename}_5_column.bed -w 2 | wc -l)
#    total=$(cat $a | wc -l)
#    results=$( echo "$bound / $total" | bc -l )
#    echo "$cluster $clusterTotal $cell $chip $investigator $bound $total $results" 
#    echo "$cluster $clusterTotal $cell $chip $investigator $bound $total $results" >> atac_chip_overlap_clusters_0.txt
#  done
#done

for a in `ls $dir0/Analysis/ATAC_seq/Union/R_output_Common_Merge0/Clusters/*15-clusters.bed`; do
  location=$(echo "$a" | rev | cut -d'/' -f2 | rev)
  cluster=$(echo "$a" | rev | cut -d'/' -f1 | rev | cut -d'.' -f3)
  clusterTotal=$(echo "$a" | rev | cut -d'/' -f1 | rev | cut -d'.' -f5 | cut -d'-' -f1)
  lines=$(wc -l $infodir/chip_seq_list.txt | cut -d' ' -f1)
  for ((i=1; i<=$lines; i++)); do
    line=$(sed -n "${i}p" < $infodir/chip_seq_list.txt)
    peakFile=$(echo "$line")
    fullName=$(echo "$peakFile" | rev | cut -d'/' -f1 | rev)
    filename=${fullName%_tag_peaks.bedgraph}
    cell=$(echo "$filename" | cut -d'_' -f1)
    chip=$(echo "$filename" | cut -d'_' -f4)
    investigator=$(echo "$filename" | cut -d'_' -f5)
    bound=$(bedtools window -a $a -b ${filename}_5_column.bed -w 2 | wc -l)
    awk -v cluster=$cluster 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"cluster "\t+"}' $a | bedtools intersect -a - -b ${filename}_5_column.bed -wb > ${filename}_${cluster}_0_merge.bed
    total=$(cat $a | wc -l)
    results=$( echo "$bound / $total" | bc -l )
    if [[ $(cat ${filename}_${cluster}_0_merge.bed | wc -l) > 0 ]]; then 
    average=$(awk '{sum+=$11} END { print sum/NR}' ${filename}_${cluster}_0_merge.bed)
    else
    average=0
    fi
    echo "$cluster $clusterTotal $cell $chip $investigator $bound $total $results $average"
    echo "$cluster $clusterTotal $cell $chip $investigator $bound $total $results $average" >> atac_chip_overlap_clusters_0_merge.txt
  done
done

#sed ' ' atac_chip_overlap.txt 
#echo -e "$directory/t$filename $(cat ${filename}_5_column.bed | wc -l) (bedtools window -a $dir1/$filename/${filename}.bed -b ${filename}_5_column.bed -w 2 | wc -l) ETS1_DP $(bedtools window -a $dir1/$filename/${filename}.bed -b $dir2/other/GSE56393_ETS1_DP_mm10.bedgraph -w 2 | wc -l) TCF1 $(bedtools window -a $dir1/$filename/${filename}.bed -b $dir2/macs/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph -w 2 | wc -l)" >> $dir1/tcf1_ets1.txt

