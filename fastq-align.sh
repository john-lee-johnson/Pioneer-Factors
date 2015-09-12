#!/bin/bash
dir=/mnt/data1/John/pioneer/data

function star_generate {
STAR --runMode genomeGenerate --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38_81 --genomeFastaFiles /mnt/data0/John/Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile /mnt/data0/John/Mus_musculus.GRCm38.81.gtf
}

function star_chip {
STAR --runMode alignReads --alignIntronMax 1 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38_81 --readFilesIn ${fastq}.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
mv Aligned.sortedByCoord.out.bam ${filename}.bam
}

function star_rna {
STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38_81 --readFilesIn ${fastq}.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
mv Aligned.sortedByCoord.out.bam ${filename}.bam
}

while read line; do
srr=$(echo $line | cut -d' ' -f1)
seq=$(echo $line | cut -d' ' -f2 | cut -d'_' -f1)
cell=$(echo $line | cut -d' ' -f2 | cut -d'_' -f2)
if [[ "$seq" = ATAC ]]; then
  cd $dir/atac_seq
  mkdir -p $cell
  cd $cell
  #fastq-dump $srr
  echo $srr $seq $cell
fi
cd $dir
if [[ "$seq" = H3* ]]
then
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
  cd $dir/rna_seq
  mkdir -p $cell
  cd $cell
  #fastq-dump $srr
  echo $srr $seq $cell
fi
cd $dir
done < files/srr_files.txt

cd /mnt/data1/John/Pioneer-Factors
while read line; do
seq=$(echo $line | cut -d' ' -f1)
cell=$(echo $line | cut -d' ' -f2)
if [[ "$seq" = ATAC ]]; then
  cd $dir/atac_seq/$cell
  #cat SRR*.fastq > ${cell}_${seq}.fastq
  fastq=${cell}_${seq}.fastq
  filename=${cell}_${seq}
  #pigz ${fastq}
  #star_chip
  rm SRR*.fastq
  echo $seq $cell
fi
cd $dir
if [[ "$seq" = H3* ]]
then
  cd $dir/chip_seq/$cell/$seq
  #cat SRR*.fastq > ${cell}_${seq}.fastq
  fastq=${cell}_${seq}.fastq
  filename=${cell}_${seq}
  #pigz ${fastq}
  star_chip
  rm SRR*.fastq
  echo $seq $cell
fi
cd $dir
if [[ "$seq" = RNA ]]
then
  cd $dir/rna_seq/$cell
  #cat SRR*.fastq > ${cell}_${seq}.fastq
  fastq=${cell}_${seq}.fastq
  filename=${cell}_${seq}
  #pigz ${fastq}
  star_rna
  rm SRR*.fastq
  echo $seq $cell
fi
cd $dir
done < files/sample_files.txt