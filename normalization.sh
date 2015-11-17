#!/bin/bash

#Getting total number of tags for ATAC-seq peaks used for normalization later
rm -f $infodir/atac_read_count.txt
cd ${dir0}/Data/ATAC_seq/
for i in `ls */Log.final.out`; do
cell=$(echo $i | cut -d'/' -f1)
echo -e $cell $(awk 'BEGIN { FS = "|\t" } {if (NR == 9) {printf "%.5f\n", $2/1000000}}' $i)
done >> $infodir/atac_read_count.txt
	
#Getting total number of tags for ChIP-seq peaks used for normalization later
rm -f $infodir/chip_read_count.txt
cd ${dir0}/Data/ChIP_seq/
for i in `ls */*/Log.final.out`; do
cell=$(echo $i | cut -d'/' -f1)
mark=$(echo $i | cut -d'/' -f2)
echo -e $cell $mark $(awk 'BEGIN { FS = "|\t" } {if (NR == 9) {printf "%.5f\n", $2/1000000}}' $i)
done >> $infodir/chip_read_count.txt

