#!/bin/bash

#---------------------Setting Directories-------------------------------------------------
#maindir=$1
maindir=/mnt/data1/John/Pioneer_Factors

#---------------------HOMER MOTIF ANALYSIS------------------------------------------------
cp ${maindir}/Analysis/ChIP_seq/homer/*_5_column.bed ${maindir}/Analysis/ChIP_seq/motif
cd ${maindir}/Analysis/ChIP_seq/motif

#batchFindMotifsGenome.pl mm9 -size given -p 35 -f *_5_column.bed	##Does motif analysis at unique ATAC seq peaks
mkdir ${maindir}/Analysis/ChIP_seq/motif/MPP
cp ${maindir}/Analysis/ChIP_seq/R_output/*MPP*.bedgraph ${maindir}/Analysis/ChIP_seq/motif/MPP
cd ${maindir}/Analysis/ChIP_seq/motif/MPP
#batchFindMotifsGenome.pl mm9 -size given -p 35 -f MPP*	##Does motif analysis at unique ATAC seq peaks
batchFindMotifsGenome.pl mm9 -size given -p 35 -f *_MPP*	##Does motif analysis at unique ATAC seq peaks