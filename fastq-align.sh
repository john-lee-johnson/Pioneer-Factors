#!/bin/bash
dir=/mnt/data1/John/pioneer/data
while read line; do ##Lists the ATAC seq bam files to do peaking calling
srr=$(echo $line | cut -d' ' -f1)
seq=$(echo $line | cut -d' ' -f2 | cut -d'_' -f1)
cell=$(echo $line | cut -d' ' -f2 | cut -d'_' -f2)
if [[ "$seq" = ATAC ]]
then
  echo $srr $seq $cell
fi
if [[ "$seq" = H3* ]]
then
  echo $srr $seq $cell
fi
if [[ "$seq" = RNA ]]
then
  echo $srr $seq $cell
fi
#STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/genome_index --readFilesIn ${fastq}.fastq --outSAMtype BAM SortedByCoordinate
#Aligned.sortedByCoord.out.bam $filename.bam
done < files/srr_files.txt
