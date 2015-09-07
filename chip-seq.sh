#!/bin/bash
datadir=/mnt/data0/ncbi
cd $datadir
mkdir -p Tcf1_ChIP_seq
cd Tcf1_ChIP_seq
mkdir -p CD8_Tcf1
mkdir -p CD8_IgG
mkdir -p Thy_Tcf1
mkdir -p Thy_Input
mkdir -p DP_Tcf1
##Defining STAR alignment function
star (){
  STAR --runMode alignReads --runThreadN 40 --alignIntronMax 1 --outFilterMultimapNmax 1 --genomeDir /mnt/data0/John/genome_index --readFilesIn $fastq --outSAMtype BAM SortedByCoordinate --outWigType wiggle
}
cd CD8_Tcf1
fastq-dump SRR1024054
fastq=SRR1024054
star

cd CD8_IgG
fastq-dump SRR1024055
fastq=SRR1024055
star

