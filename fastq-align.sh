#!/bin/bash
##This script will read the SRR files from srr_files.txt, download and convert them to fastq, align them, and output to BAM files


dir0=/mnt/data1/John/pioneer
filedir=/mnt/data1/John/Pioneer-Factors/files

#Makes the directories
mkdir -p $dir0/data/atac_seq/bam
mkdir -p $dir0/data/chip_seq/bam
mkdir -p $dir0/data/rna_seq/bam
mkdir -p $dir0/data/chip_seq/tcf1/bam


##---------------------SETTING FUNCTIONS---------------------------------------------------
#STAR generates genome index
function star_generate {
STAR --runMode genomeGenerate --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --genomeFastaFiles /mnt/data0/John/GRCm38p4_mm10_pa_only.fa --sjdbGTFfile /mnt/data0/John/gencode.vM6.annotation.gtf
}

#STAR alignment for ChIP-seq and ATAC-seq
function star_chip {
STAR --runMode alignReads --alignIntronMax 1 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${fastq}.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
mv Aligned.sortedByCoord.out.bam ${filename}.bam
}

#STAR alignment for RNA-seq 2-pass
function star_rna {
STAR --runMode alignReads --twopassMode Basic --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${fastq}.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
mv Aligned.sortedByCoord.out.bam ${filename}.bam
}

#STAR alignment for color space ChIP-seq
function star_color {
STAR --runMode alignReads --alignIntronMax 1 --outFilterMultimapNmax 2 --runThreadN 40 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${fastq}.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
mv Aligned.sortedByCoord.out.bam ${filename}.bam
}

#Read in a file and download the SRR file and convert to fastq
while read line; do
srr=$(echo $line | cut -d' ' -f1)
seq=$(echo $line | cut -d' ' -f2 | cut -d'_' -f1)
cell=$(echo $line | cut -d' ' -f2 | cut -d'_' -f2)
if [[ "$seq" = ATAC ]]; then
  echo $srr $seq $cell
  cd $dir0/data/atac_seq ; mkdir -p $cell ; cd $cell
  #fastq-dump $srr
fi
if [[ "$seq" = H3* ]]
then
  echo $srr $seq $cell
  cd $dir0/data/chip_seq ; mkdir -p $cell ; cd $cell ; mkdir -p $seq
  cd $seq
  #fastq-dump $srr
fi
if [[ "$seq" = RNA ]]
then
  echo $srr $seq $cell
  cd $dir0/data/rna_seq ; mkdir -p $cell ; cd $cell
  #fastq-dump $srr
fi
done < $filedir/srr_files.txt

#Read in a file and download the SRR file and convert to fastq for Tcf1
MACSpvalue=1e-7
while read line; do
  srr=$(echo $line | cut -d' ' -f1)
  sample=$(echo $line | cut -d' ' -f2 | cut -d'_' -f1)
  filename=$(echo $line | cut -d' ' -f2)
  fastq=${filename}.fastq
  cell=$(echo $line | cut -d' ' -f2 | cut -d'_' -f2)
  color=$(echo $line | cut -d' ' -f3)
  treat=$(echo $line | cut -d' ' -f4)
  cd $dir0/data/chip_seq/tcf1 ; mkdir -p $cell ; cd $cell ; mkdir -p $sample ; cd $sample
  if [ $color = "Color" ]; then
  #fastq-dump -B $srr && mv ${srr}.fastq ${filename}.fastq && pigz ${filename}.fastq
  star_color && mv ${filename}.bam $dir0/data/chip_seq/tcf1/bam
  readLength=$(awk 'BEGIN { FS = "|\t" } {if (NR == 7) {print $2}}' Log.final.out)
    if [ $treat = "Treat" ]; then
    echo "${dir0}/data/chip_seq/tcf1/bam/${filename}.bam ${filename}_macs ${MACSpvalue} $readLength $control"
    else
    control=${dir0}/data/chip_seq/tcf1/bam/${filename}.bam
    fi
  else
  #fastq-dump $srr && mv ${srr}.fastq ${filename}.fastq && pigz ${filename}.fastq
  star_chip && mv ${filename}.bam $dir0/data/chip_seq/tcf1/bam
  readLength=$(awk 'BEGIN { FS = "|\t" } {if (NR == 7) {print $2}}' Log.final.out)
    if [ $treat = "Treat" ]; then
    echo "${dir0}/data/chip_seq/tcf1/bam/${filename}.bam ${filename}_macs ${MACSpvalue} $readLength $control"
    else
    control=${dir0}/data/chip_seq/tcf1/bam/${filename}.bam
    fi
  fi
done < $filedir/tcf1.txt > $filedir/macs_parallel_tcf1.txt

#Combine fastq files, compress to .gz, and align
cd /mnt/data1/John/Pioneer-Factors
while read line; do
  seq=$(echo $line | cut -d' ' -f1)
  cell=$(echo $line | cut -d' ' -f2)
  fastq=${cell}_${seq}.fastq
  filename=${cell}_${seq}
if [[ "$seq" = ATAC ]]; then
  cd $dir0/data/atac_seq/$cell
  #cat SRR*.fastq > ${filename}.fastq && pigz ${fastq} && rm SRR*.fastq
  star_chip && mv ${filename}.bam $dir0/data/atac_seq/bam
  readLength=$(awk 'BEGIN { FS = "|\t" } {if (NR == 7) {print $2}}' Log.final.out)
  echo "${dir0}/data/atac_seq/bam/${filename}.bam ${filename}_macs ${MACSpvalue} $readLength"
fi
cd $dir0
if [[ "$seq" = H3* ]]
then
  cd $dir0/data/chip_seq/$cell/$seq
  #cat SRR*.fastq > ${cell}_${seq}.fastq && pigz ${fastq} && rm SRR*.fastq
  star_chip && mv ${cell}_${seq}.bam $dir0/data/chip_seq/bam
fi
cd $dir0
if [[ "$seq" = RNA ]]
then
  cd $dir0/data/rna_seq/$cell
  #cat SRR*.fastq > ${cell}_${seq}.fastq && pigz ${fastq} && rm SRR*.fastq
  star_rna && mv ${cell}_${seq}.bam $dir0/data/rna_seq/bam
fi
cd $dir0
done < $filedir/sample_files.txt > $filedir/macs_parallel_atac.txt
