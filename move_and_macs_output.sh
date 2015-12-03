#!/bin/bash

rm -f $paralleldir/macs_parallel_atac_srr*.txt
rm -f $paralleldir/macs_parallel_chip_srr*.txt
rm -f $paralleldir/homer_key_file.txt
rm -f $paralleldir/macs_parallel_atac.txt
rm -f $paralleldir/macs_parallel_chip.txt

##---------------------MOVE DATA FILES AMIT-----------------------------------------------
#Combine fastq files, compress to .gz, and align
lines=$(wc -l $infodir/sample_files_amit.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_amit.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
seq=$(echo $line | cut -d':' -f2 | cut -d';' -f3 | xargs)
cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq/$cell
  filename=${cell}_${seq}_Amit_mm10_${gsm}
  samtools index ${filename}.bam
  cp ${filename}.bam $datadir/ATAC_seq/bam
  cp ${filename}.bam.bai $datadir/ATAC_seq/bam
  cp ${filename}.fastq.gz $datadir/ATAC_seq/fastq
  echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
fi
if [[ "$seq" = ChIP* ]]; then
  mark=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs | cut -d'_' -f1 | xargs)
  cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs | cut -d'_' -f2 | xargs)
  filename=${cell}_${seq}_${mark}_Amit_mm10_${gsm}
  cd $dir0/Data/ChIP_seq/$cell/$mark
  samtools index ${filename}.bam
  cp ${filename}.bam $datadir/ChIP_seq/bam
  cp ${filename}.bam.bai $datadir/ChIP_seq/bam
  samtools index $datadir/ChIP_seq/bam/${filename}.bam
  cp ${filename}.fastq.gz $datadir/ChIP_seq/fastq
  echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq/$cell
  filename=${cell}_${seq}_Amit_mm10_${gsm}
  samtools index ${filename}.bam
  cp ${filename}.bam $datadir/RNA_seq/bam
  cp ${filename}.bam.bai $datadir/RNA_seq/bam
  samtools index $datadir/RNA_seq/bam/${filename}.bam
  cp ${filename}.fastq.gz $datadir/RNA_seq/fastq
fi
done

##--------------MACS OUTPUT COMMANDS AMIT SEPARATE SRR FILES------------------------------
#Align fastq to bam for Amit data
while read line; do
seq=$(echo $line | cut -d' ' -f1)
if [[ "$seq" = ATAC* ]]; then
  cell=$(echo $line | cut -d' ' -f2 | xargs)
  cd $dir0/Data/ATAC_seq/$cell
  mkdir -p ${dir0}/Analysis/ATAC_seq/srr_macs/$cell
  srr_line=${line##"$seq $cell "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo `pwd`"/${srr}.bam ${srr}_macs" >> $paralleldir/macs_parallel_atac_srr_${cell}.txt
  done
fi
if [[ "$seq" = ChIP* ]]; then
  mark=$(echo $line | cut -d' ' -f2 | xargs)
  cell=$(echo $line | cut -d' ' -f3 | xargs)
  cd $dir0/Data/ChIP_seq/$cell/$mark
  mkdir -p ${dir0}/Analysis/ChIP_seq/srr_macs/$cell/$mark
  srr_line=${line##"$seq $mark $cell "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo `pwd`"/${srr}.bam ${srr}_macs" >> $paralleldir/macs_parallel_chip_srr_${cell}_${mark}.txt
  done
fi
if [[ "$seq" = RNA* ]]
then
  cell=$(echo $line | cut -d' ' -f2)
  cd $dir0/Data/RNA_seq/$cell
  srr_line=${line##"$seq $cell "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
  echo $srr
    #mark_duplicates #consensus is not to remove duplicates from RNA-seq data
  done
fi
done < $infodir/srr_files_amit.txt

cd $dir0/Analysis/Homer/Tag_Directories/Amit
batchMakeTagDirectory.pl $paralleldir/homer_key_file.txt -cpu 40 -genome mm10 -format sam

##--------------MACS OUTPUT COMMANDS AMIT COMBINED BAM FILE-------------------------------
#Read in a file and output commands to a file for MACS peak calling for Amit Data
lines=$(wc -l $infodir/sample_files_amit.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_amit.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
seq=$(echo $line | cut -d':' -f2 | cut -d';' -f3 | xargs)
cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq/$cell
  filename=${cell}_${seq}_Amit_mm10_${gsm}
  fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/Amit/$filename/tagInfo.txt" | cut -d"=" -f2 | xargs)
  echo "$datadir/ATAC_seq/bam/${filename}.bam ${filename}_macs" "$fragLength">> $paralleldir/macs_parallel_atac.txt
fi
if [[ "$seq" = ChIP* ]]; then
  mark=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | cut -d'_' -f1 | xargs)
  cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | cut -d'_' -f2 | xargs)
  filename=${cell}_${seq}_${mark}_Amit_mm10_${gsm}
  fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/Amit/$filename/tagInfo.txt" | cut -d"=" -f2 | xargs)
  echo "$datadir/ChIP_seq/bam/${filename}.bam ${filename}_macs" "$fragLength">> $paralleldir/macs_parallel_chip.txt
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq/$cell
fi
done