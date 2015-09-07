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
  STAR --runMode alignReads --runThreadN 40 --alignIntronMax 1 --outFilterMultimapNmax 1 --genomeDir /mnt/data0/John/GRCh38.81.mm --readFilesIn $fastq --outSAMtype BAM SortedByCoordinate
}
cd $datadir/Tcf1_ChIP_seq/CD8_Tcf1
fastq-dump SRR1024054
fastq=SRR1024054.fastq
star

cd $datadir/Tcf1_ChIP_seq/CD8_IgG
fastq-dump SRR1024055
fastq=SRR1024055.fastq
star

cd $datadir/Tcf1_ChIP_seq/Thy_Tcf1
fastq-dump SRR846897
fastq=SRR846897.fastq
star

cd $datadir/Tcf1_ChIP_seq/Thy_Input
fastq-dump SRR846898
fastq=SRR846898.fastq
star

cd $datadir/Tcf1_ChIP_seq/DP_Tcf1
fastq-dump SRR1685994
fastq=SRR1685994.fastq
star
