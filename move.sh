#!/bin/bash
cd /mnt/data1/John/pioneer/data/chip_seq/bam
for i in `ls *bam`; do
filename=$(echo $i | cut -d'.' -f1)
#echo $filename
cp $i /mnt/data1/VahediLab/IdoAmit_Science_Data/Data/ChIP-seq/bam_JJ/${filename}_mm10.bam
done

cd /mnt/data1/John/pioneer/data/chip_seq/
for i in `ls */*/*fastq.gz`; do
filename=$(echo $i | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
echo $filename
#cp $i /mnt/data1/VahediLab/IdoAmit_Science_Data/Data/ChIP-seq/fastq_JJ/${filename}_mm10.fastq.gz
done
