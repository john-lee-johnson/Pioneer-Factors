#!/bin/bash

#---------------------Setting Directories-------------------------------------------------
#maindir=$1
maindir=/mnt/data1/John/pioneer

#---------------------HOMER MOTIF ANALYSIS------------------------------------------------
cd ${maindir}/analysis/chip_seq/tcf1/macs/motif
for i in `ls *macs_peaks*`; do
changeNewLine.pl $i
done
batchFindMotifsGenome.pl mm10 -size given -p 35 -f *macs_peaks* ##Does motif analysis at Tcf1 peaks
