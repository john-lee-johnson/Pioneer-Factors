#!/bin/bash
while read line; do ##Lists the ATAC seq bam files to do peaking calling
srr=$(echo $line | cut -d' ' -f1)
filename=$(echo $line | cut -d' ' -f2)
echo $srr
echo $filename
#STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/genome_index --readFilesIn ${fastq}.fastq --outSAMtype BAM SortedByCoordinate
Aligned.sortedByCoord.out.bam $filename.bam
done < fastq_files.txt
