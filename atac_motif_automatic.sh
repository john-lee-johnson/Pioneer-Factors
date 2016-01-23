#!/bin/bash
dir0=/mnt/data1/John/Pioneer_Factors
cd $dir0/Analysis/ATAC_seq/Union/R_output_Common_0.5
for a in `ls -d */`; do
  cd $a
  filename=$(echo $a | cut -d'/' -f1)
  cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed .
  cp $dir0/Analysis/ChIP_seq/Transcription_Factor/Motif/Tcf1/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed .
  mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
  rm -f Tcf1.txt
  bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b ${filename}.bed -wa > ${filename}_peaks.bed
  bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b ${filename}_peaks.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > ${filename}_peaks_subtract.bed
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${filename}_peaks_subtract.bed > ${filename}_peaks_subtract_5_column.bed ##Adds a unique peak name and strand info
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${filename}.bed > ${filename}_peaks_5_column_merge.bed
  findMotifsGenome.pl ${filename}_peaks_5_column_merge.bed mm10 Motif.$filename -size 200 -p 35 -bg ${filename}_peaks_subtract_5_column.bed -keepFiles	##Does motif analysis at unique ATAC seq peaks
  bedtools intersect -a ${filename}_peaks_5_column_merge.bed -b Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed > Tcf1_peaks.bed
  tcf1=$(wc -l Tcf1_peaks.bed | cut -d" " -f1)
  peaks=$(wc -l ${filename}_peaks_5_column_merge.bed | cut -d" " -f1)
  results=$( echo "$tcf1 / $peaks" | bc -l )
  if [[ $(echo "if (${results} > 0.6) 1 else 0" | bc) -eq 1 ]]; then
    #findMotifsGenome.pl ${filename}_peaks_5_column_merge.bed mm10 Motif.$filename -bg ${filename}_peaks_subtract_5_column.bed -nomotif -dumpFasta -size 200 -p 35 -keepFiles	##Does motif analysis at unique ATAC seq peaks
    #findMotifsGenome.pl ${filename}_peaks_5_column_merge.bed mm10 Motif.$filename -size 200 -p 35 -bg ${filename}_peaks_subtract_5_column.bed -keepFiles	##Does motif analysis at unique ATAC seq peaks
    count=`head Motif.$filename/knownResults.txt -n 5 | grep -i Tcf | grep -v CTCF | grep -v Tcf12 | wc -l`
    echo "Greater than 60 percent! FC < 1.5 $results Motif.$filename"
    if [ $count -ne 0 ]
    then
     line=$(head Motif.$filename/knownResults.txt -n 5 | grep -i Tcf | grep -v CTCF | grep -v Tcf12)
      echo "$line" >> Tcf1.txt
      echo "$line FC < 1.5 $results Motif.$filename" | mail -s "Tcf1 found FC" johnlee.johnson@gmail.com
      unset count
      echo "MATCH FOUND! Motif.$num"
    else
      echo "No match with the pattern"
    fi
    else
     echo "Less than 60 percent!"
  fi
  cd $dir0/Analysis/ATAC_seq/Union/R_output_Common_0.5
done

dir0=/mnt/data1/John/Pioneer_Factors
cd $dir0/Analysis/ATAC_seq/Union/R_output_Common_1.5
for a in `ls -d */`; do
  cd $a
  filename=$(echo $a | cut -d'/' -f1)
  cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed .
  cp $dir0/Analysis/ChIP_seq/Transcription_Factor/Motif/Tcf1/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed .
  mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
  rm -f Tcf1.txt
  bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b ${filename}.bed -wa > ${filename}_peaks.bed
  bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b ${filename}_peaks.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > ${filename}_peaks_subtract.bed
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${filename}_peaks_subtract.bed > ${filename}_peaks_subtract_5_column.bed ##Adds a unique peak name and strand info
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${filename}.bed > ${filename}_peaks_5_column_merge.bed
  findMotifsGenome.pl ${filename}_peaks_5_column_merge.bed mm10 Motif.$filename -size 200 -p 35 -bg ${filename}_peaks_subtract_5_column.bed -keepFiles
  bedtools intersect -a ${filename}.bed -b Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed > Tcf1_peaks.bed
  tcf1=$(wc -l Tcf1_peaks.bed | cut -d" " -f1)
  peaks=$(wc -l ${filename}_peaks_5_column_merge.bed | cut -d" " -f1)
  results=$( echo "$tcf1 / $peaks" | bc -l )
  if [[ $(echo "if (${results} > 0.59) 1 else 0" | bc) -eq 1 ]]; then
    #findMotifsGenome.pl ${filename}_peaks_5_column_merge.bed mm10 Motif.$filename -bg ${filename}_peaks_subtract_5_column.bed -nomotif -dumpFasta -size 200 -p 35 -keepFiles	##Does motif analysis at unique ATAC seq peaks
    #findMotifsGenome.pl ${filename}_peaks_5_column_merge.bed mm10 Motif.$filename -size 200 -nomotif -p 35 -bg ${filename}_peaks_subtract_5_column.bed -keepFiles
    count=`head Motif.$filename/knownResults.txt -n 5 | grep -i Tcf | grep -v CTCF | grep -v Tcf12 | wc -l`
    echo "Greater than 60 percent! FC > 1.5 $results Motif.$filename"
    if [ $count -ne 0 ]
    then
      line=$(head Motif.$filename/knownResults.txt -n 5 | grep -i Tcf | grep -v CTCF | grep -v Tcf12)
      echo "$line" >> Tcf1.txt
      echo "$line  FC > 1.5 $results Motif.$filename" | mail -s "Tcf1 found" johnlee.johnson@gmail.com
      unset count
      echo "MATCH FOUND! Motif.$num"
    else
      echo "No match with the pattern"
    fi
  fi
  cd $dir0/Analysis/ATAC_seq/Union/R_output_Common_1.5
done

dir0=/mnt/data1/John/Pioneer_Factors
cd $dir0/Analysis/ATAC_seq/Union/R_output_Common_2
for a in `ls -d */`; do
  cd $a
  filename=$(echo $a | cut -d'/' -f1)
  cp $dir0/Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed .
  cp $dir0/Analysis/ChIP_seq/Transcription_Factor/Motif/Tcf1/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed .
  mergeBed -i B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed > B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed
  rm -f Tcf1.txt
  bedtools intersect -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b ${filename}.bed -wa > ${filename}_peaks.bed
  bedtools subtract -A -a B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed -b ${filename}_peaks.bed | awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' - > ${filename}_peaks_subtract.bed
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${filename}_peaks_subtract.bed > ${filename}_peaks_subtract_5_column.bed ##Adds a unique peak name and strand info
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t0" "\t+"}' ${filename}.bed > ${filename}_peaks_5_column_merge.bed
  findMotifsGenome.pl ${filename}_peaks_5_column_merge.bed mm10 Motif.$filename -size 200 -p 35 -bg ${filename}_peaks_subtract_5_column.bed -keepFiles  
  bedtools intersect -a ${filename}_peaks_5_column_merge.bed -b Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed > Tcf1_peaks.bed
  tcf1=$(wc -l Tcf1_peaks.bed | cut -d" " -f1)
  peaks=$(wc -l ${filename}_peaks_5_column_merge.bed | cut -d" " -f1)
  results=$( echo "$tcf1 / $peaks" | bc -l )
  if [[ $(echo "if (${results} > 0.6) 1 else 0" | bc) -eq 1 ]]; then
    #findMotifsGenome.pl ${filename}_peaks_5_column_merge.bed mm10 Motif.$filename -bg ${filename}_peaks_subtract_5_column.bed -nomotif -dumpFasta -size 200 -p 35 -keepFiles	##Does motif analysis at unique ATAC seq peaks
    #findMotifsGenome.pl ${filename}_peaks_5_column_merge.bed mm10 Motif.$filename -size 200 -p 35 -bg ${filename}_peaks_subtract_5_column.bed -keepFiles	##Does motif analysis at unique ATAC seq peaks
    count=`head Motif.$filename/knownResults.txt -n 5 | grep -i Tcf | grep -v CTCF | grep -v Tcf12 | wc -l`
    echo "Greater than 60 percent! FC > 2 $results Motif.$filename"
    if [ $count -ne 0 ]
    then
      line=$(head Motif.$filename/knownResults.txt -n 5 | grep -i Tcf | grep -v CTCF | grep -v Tcf12)
      echo "$line" >> Tcf1.txt
      echo "$line  FC > 2 $results Motif.$filename" | mail -s "Tcf1 found" johnlee.johnson@gmail.com
      unset count
      echo "MATCH FOUND! Motif.$num"
    else
      echo "No match with the pattern"
    fi
  fi
  cd $dir0/Analysis/ATAC_seq/Union/R_output_Common_2
done