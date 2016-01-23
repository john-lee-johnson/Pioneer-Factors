#!/bin/bash
dir0=/mnt/data1/John/Pioneer_Factors
cd $dir0/Analysis/ATAC_seq/Union/R_output_Discover/Cluster_Plots/
for a in `ls -d K*`; do
  cd $a
  cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed .
  cp $dir0/Analysis/ChIP_seq/Transcription_Factor/Motif/Tcf1/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed .
  mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
  rm -f Tcf1.txt
  for i in `ls Clustered*bed`; do
    num=$(echo "$i" | cut -d'.' -f3)
    bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b $i -wa > peaks.bed
    bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b $i | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > peaks_subtract_merge.bed
    awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' peaks_subtract_merge.bed > peaks_subtract_merge_5_column.bed ##Adds a unique peak name and strand info
    awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' peaks.bed > peaks_5_column_merge.bed
    #findMotifsGenome.pl peaks_5_column_merge.bed mm10 Motif.$num -size given -nomotif -p 35 -bg peaks_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks
    bedtools intersect -a peaks_5_column_merge.bed -b Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed > Tcf1_peaks.bed
    tcf1=$(wc -l Tcf1_peaks.bed | cut -d" " -f1)
    peaks=$(wc -l peaks_5_column_merge.bed | cut -d" " -f1)
    #echo print "$tcf1"/"$peaks". | python
    results=$( echo "$tcf1 / $peaks" | bc -l )
    echo "$results Motif.$num""
    if [[ $(echo "if (${results} > 0.7) 1 else 0" | bc) -eq 1 ]]; then
      findMotifsGenome.pl peaks_5_column_merge.bed mm10 Motif.$num -size given -nomotif -p 35 -bg peaks_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks
      count=`head Motif.${num}/knownResults.txt -n 5 | grep -i Tcf | grep -v CTCF | grep -v Tcf12 | wc -l`
      echo "Greater than 70 percent!"
      if [ $count -ne 0 ]
      then
        line=$(head Motif.${num}/knownResults.txt -n 5 | grep -i Tcf | grep -v CTCF | grep -v Tcf12)
        echo "$line" >> Tcf1.txt
        echo "$line" | mail -s "Tcf1 found" johnlee.johnson@gmail.com
        unset count
        echo "MATCH FOUND! Motif.$num"
     else
        echo "No match with the pattern"
     fi
    else
     echo "Less than 60 percent!"
    fi
  done
  cd $dir0/Analysis/ATAC_seq/Union/R_output_Discover/Cluster_Plots/
done

#count=`head Motif.${num}/knownResults.txt | grep -i Tcf | grep -v CTCF | grep -v Tcf12 | wc -l`
#if [ $count -ne 0 ]
#then
#     line=$(head Motif.${num}/knownResults.txt | grep -i Tcf | grep -v CTCF | grep -v Tcf12)
#     echo "$line" >> Tcf1.txt
#     echo "$line" | mail -s "Tcf1 found" johnlee.johnson@gmail.com
#else
#     echo "no match with the pattern"
#fi
#unset count


#rm -f Tcf1.txt
#for i in `ls */knownResults.txt`; do
#count=`head $i | grep -i Tcf | wc -l`
#if [ $count -ne 0 ]
#then
#     line=$(head $i | grep -i Tcf)
#     echo "$line" >> Tcf1.txt
#     echo "$line" | mail -s "Tcf1 found" johnlee.johnson@gmail.com
#else
#     echo "no match with the pattern"
#fi
#unset count
#num=$(echo "$i" | cut -d'.' -f3)
#bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b $i -wa > peaks.bed
#bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b $i | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > peaks_subtract_merge.bed

#awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' peaks_subtract_merge.bed > peaks_subtract_merge_5_column.bed ##Adds a unique peak name and strand info
#awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' peaks.bed > peaks_5_column_merge.bed
#findMotifsGenome.pl peaks_5_column_merge.bed mm10 Motif.$num -size given -nomotif -p 35 -bg peaks_subtract_merge_5_column.bed	##Does motif analysis at unique ATAC seq peaks

#done

