#!/bin/bash
##This script will read the SRR files from srr_files.txt, download and convert them to fastq, align them, and output to BAM files

dir=/mnt/data1/John/pioneer/data

#generates genome index for STAR
#change chr names in genome index chrName.txt to chr1, chr2, chr3...
function star_generate {
STAR --runMode genomeGenerate --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38_81 --genomeFastaFiles /mnt/data0/John/Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile /mnt/data0/John/Mus_musculus.GRCm38.81.gtf
}

#STAR alignment for ChIP-seq and ATAC-seq
function star_chip {
STAR --runMode alignReads --alignIntronMax 1 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38_81 --readFilesIn ${fastq}.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
mv Aligned.sortedByCoord.out.bam ${filename}.bam
}

#STAR alignment for RNA-seq 2-pass
function star_rna {
STAR --runMode alignReads --twopassMode Basic --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38_81 --readFilesIn ${fastq}.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
mv Aligned.sortedByCoord.out.bam ${filename}.bam
}

#Read in a file and download the SRR file and convert to fastq
while read line; do
srr=$(echo $line | cut -d' ' -f1)
seq=$(echo $line | cut -d' ' -f2 | cut -d'_' -f1)
cell=$(echo $line | cut -d' ' -f2 | cut -d'_' -f2)
if [[ "$seq" = ATAC ]]; then
  echo $srr $seq $cell
  cd $dir/atac_seq
  mkdir -p $cell
  cd $cell
  #fastq-dump $srr
fi
cd $dir
if [[ "$seq" = H3* ]]
then
  echo $srr $seq $cell
  cd $dir/chip_seq
  mkdir -p $cell
  cd $cell
  mkdir -p $seq
  cd $seq
  #fastq-dump $srr
fi
cd $dir
if [[ "$seq" = RNA ]]
then
  echo $srr $seq $cell
  cd $dir/rna_seq
  mkdir -p $cell
  cd $cell
  #fastq-dump $srr
fi
cd $dir
done < files/srr_files.txt

mkdir -p $dir/rna_seq/bam
mkdir -p $dir/chip_seq/bam
mkdir -p $dir/atac_seq/bam

#Combine fastq files, compress to .gz, and align
cd /mnt/data1/John/Pioneer-Factors
while read line; do
seq=$(echo $line | cut -d' ' -f1)
cell=$(echo $line | cut -d' ' -f2)
if [[ "$seq" = ATAC ]]; then
  echo $seq $cell
  cd $dir/atac_seq/$cell
  #cat SRR*.fastq > ${cell}_${seq}.fastq
  fastq=${cell}_${seq}.fastq
  filename=${cell}_${seq}
  #pigz ${fastq}
  #star_chip
  #rm SRR*.fastq
  #mv ${cell}_${seq}.bam $dir/atac_seq/bam
fi
cd $dir
if [[ "$seq" = H3* ]]
then
  echo $seq $cell
  cd $dir/chip_seq/$cell/$seq
  #cat SRR*.fastq > ${cell}_${seq}.fastq
  fastq=${cell}_${seq}.fastq
  filename=${cell}_${seq}
  #pigz ${fastq}
  #star_chip
  #rm SRR*.fastq
  #mv ${cell}_${seq}.bam $dir/atac_seq/bam
fi
cd $dir
if [[ "$seq" = RNA ]]
then
  echo $seq $cell
  cd $dir/rna_seq/$cell
  #cat SRR*.fastq > ${cell}_${seq}.fastq
  fastq=${cell}_${seq}.fastq
  filename=${cell}_${seq}
  #pigz ${fastq}
  star_rna
  #rm SRR*.fastq
  mv ${cell}_${seq}.bam $dir/rna_seq/bam
fi
cd $dir
done < files/sample_files.txt