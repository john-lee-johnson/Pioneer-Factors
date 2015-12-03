#!/bin/bash

#Getting total number of tags for ATAC-seq peaks used for normalization later
rm -f $infodir/atac_read_count.txt
cd $datadir/ATAC_seq/bam
for i in `ls *.bam`; do
cell=$(echo $i | cut -d'_' -f1)
seq=$(echo $i | cut -d'_' -f2)
pi=$(echo $i | cut -d'_' -f4)
counts=$(samtools view -c $i)
echo -e "$cell" "$seq" "$pi" $(bc <<< "scale = 3; $counts/1000000") >> $infodir/atac_read_count.txt
done
	
#Getting total number of tags for ChIP-seq peaks used for normalization later
rm -f $infodir/chip_read_count.txt
cd $datadir/ChIP_seq/bam/
for i in `ls *.bam`; do
cell=$(echo $i | cut -d'/' -f1)
mark=$(echo $i | cut -d'/' -f4)
pi=$(echo $i | cut -d'_' -f5)
counts=$(samtools view -c $i)
echo -e "$cell" "$mark" "$pi" $(bc <<< "scale = 3; $counts/1000000") >> $infodir/chip_read_count.txt
done

