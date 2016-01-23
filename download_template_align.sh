#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory

##---------------------SETTING FUNCTIONS FOR ALIGNMENT-------------------------------------
#Alignment functions should ultimately generate the aligned file in the following format:
#SRR123456.bam

#STAR generates genome index mm10
function star_generate_mm10 {
STAR --runMode genomeGenerate --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --genomeFastaFiles /mnt/data0/John/GRCm38p4_mm10_pa_only.fa --sjdbGTFfile /mnt/data0/John/gencode.vM6.annotation.gtf
}
#STAR generates genome index hg19
function star_generate_hg19 {
STAR --runMode genomeGenerate --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --genomeFastaFiles /mnt/data0/John/GRCh37.p13.genome.fa --sjdbGTFfile /mnt/data0/John/gencode.v19.annotation.gtf
}
#STAR alignment for ATAC-seq
function star_atac {
echo "STAR VERSION: `STAR --version`" >> $logdir/download_log_${des}.txt
if [[ "$genome" = mm* ]]; then
  if [[ "$type" = Pair* ]]; then
    echo "STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_1.fastq ${srr}_2.fastq --outSAMtype BAM SortedByCoordinate" >> $logdir/download_log_${des}.txt
    STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn${srr}_1.fastq ${srr}_2.fastq --outSAMtype BAM SortedByCoordinate 
  elif [[ "$type" = Single* ]]; then
    echo "STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate" >> $logdir/download_log_${des}.txt
    STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate   
  fi
elif [[ "$genome" = hg* ]]; then
  if [[ "$type" = Pair* ]]; then
    echo "STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}_1.fastq ${srr}_2.fastq --outSAMtype BAM SortedByCoordinate" >> $logdir/download_log_${des}.txt
    STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}_1.fastq ${srr}_2.fastq --outSAMtype BAM SortedByCoordinate 
  elif [[ "$type" = Single* ]]; then
    echo "STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate" >> $logdir/download_log_${des}.txt
    STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate 
  fi
fi
mv Aligned.sortedByCoord.out.bam ${srr}.bam #Keeps the bam file that STAR generates
mv Log.final.out ${srr}_Log.final.out #Keeps a copy of the aligner log
}
#ABI alignment using BOWTIE
function abi-align {
echo "BOWTIE VERSION: `bowtie --version`" >> $logdir/download_log_${des}.txt
if [[ "$genome" = mm* ]]; then
  echo "bowtie -p 38 -S -C /mnt/data0/John/bowtie_mm10_colorspace/mm10_colorspace -f ${srr}_F3.csfasta -Q ${srr}_F3_QV.qual ${srr}.sam" >> $logdir/download_log_${des}.txt
  bowtie -p 38 -S -C /mnt/data0/John/bowtie_mm10_colorspace/mm10_colorspace -f ${srr}_F3.csfasta -Q ${srr}_F3_QV.qual ${srr}.sam
elif [[ "$genome" = hg* ]]; then
  echo "bowtie -p 38 -S -C /mnt/data0/John/bowtie_hg19_colorspace/hg19_c -f ${srr}_F3.csfasta -Q ${srr}_F3_QV.qual ${srr}.sam" >> $logdir/download_log_${des}.txt
  bowtie -p 38 -S -C /mnt/data0/John/bowtie_hg19_colorspace/hg19_c -f ${srr}_F3.csfasta -Q ${srr}_F3_QV.qual ${srr}.sam
fi
mkdir -p tmp
echo "SAMTOOLS VERSION: `samtools --version`" >> $logdir/download_log_${des}.txt
samtools view -b -h ${srr}.sam > ${srr}_all_reads.bam #Convert sam to bam to save space
rm -f ${srr}.sam
echo "samtools view -b -h -F 4 ${srr}.sam > ${srr}_all_reads.bam #Keep only aligned reads" >> $logdir/download_log_${des}.txt
samtools view -b -h -F 4 ${srr}_all_reads.bam > ${srr}.bam #Keep only aligned reads
mv ${srr}.bam ${srr}_sort.bam
echo "------------SORTING BAM FILE-----------------" >> $logdir/download_log_${des}.txt #Sort bam file according to coordinates
echo "PICARD VERSION: `java -jar $PICARD SortSam --version`" >> $logdir/download_log_${des}.txt
echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD SortSam INPUT=${srr}_sort.bam OUTPUT=${srr}.bam SORT_ORDER=coordinate" >> $logdir/download_log_${des}.txt
java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD SortSam INPUT=${srr}_sort.bam OUTPUT=${srr}.bam SORT_ORDER=coordinate
#cp ${srr}_F3.csfasta $datadir/$seq/reads
#cp ${srr}_F3_QV.qual $datadir/$seq/reads
}

#STAR alignment for ChIP-seq
function star_chip {
echo "STAR VERSION: `STAR --version`" >> $logdir/download_log_${des}.txt
if [[ "$genome" = mm* ]]; then
  if [[ "$type" = Pair* ]]; then
    echo "STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_1.fastq ${srr}_2.fastq --outSAMtype BAM SortedByCoordinate" >> $logdir/download_log_${des}.txt
    STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_1.fastq ${srr}_2.fastq --outSAMtype BAM SortedByCoordinate 
  elif [[ "$type" = Single* ]]; then
    echo "STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate" >> $logdir/download_log_${des}.txt
    STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}.fastq  --outSAMtype BAM SortedByCoordinate 
  fi
elif [[ "$genome" = hg* ]]; then
  if [[ "$type" = Pair* ]]; then
    echo "STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}_1.fastq ${srr}_2.fastq --outSAMtype BAM SortedByCoordinate" >> $logdir/download_log_${des}.txt
    STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}_1.fastq ${srr}_2.fastq --outSAMtype BAM SortedByCoordinate 
  elif [[ "$type" = Single* ]]; then
    echo "STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate" >> $logdir/download_log_${des}.txt
    STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate 
  fi
fi
}

#STAR alignment for RNA-seq 2-pass
function star_rna {
echo "STAR VERSION: `STAR --version`" >> $logdir/download_log_${des}.txt
if [[ "$genome" = mm* ]]; then
  if [[ "$type" = Pair* ]]; then
    echo "STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_1.fastq ${srr}_2.fastq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30" >> $logdir/download_log_${des}.txt
    STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_1.fastq ${srr}_2.fastq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30
  elif [[ "$type" = Single* ]]; then
    echo "STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30" >> $logdir/download_log_${des}.txt
    STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30  
  fi
elif [[ "$genome" = hg* ]]; then
  if [[ "$type" = Pair* ]]; then
    echo "STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}_1.fastq ${srr}_2.fastq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30" >> $logdir/download_log_${des}.txt
    STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}_1.fastq ${srr}_2.fastq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30
  elif [[ "$type" = Single* ]]; then
    echo "STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30" >> $logdir/download_log_${des}.txt
    STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30  
  fi
fi
mv Aligned.sortedByCoord.out.bam ${srr}.bam #Keeps the bam file that STAR generates
mv Log.final.out ${srr}_Log.final.out #Keeps a copy of the aligner log
}

#STAR remove genome from memory
function star_remove {
STAR --genomeLoad Remove --genomeDir /mnt/data0/John/genome_GRCm38p4_M6
}

###---------------------ALIGNMENT---------------------------------
#Read in SRR files and align
echo ""; echo "#--------------ALIGNMENT--------------------"
echo "#--------------ALIGNMENT-------------------------" >> $logdir/download_log_${des}.txt
while read line; do
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
investigator=$(echo $line | cut -d',' -f5)
sequencer=$(echo $line | cut -d',' -f6)
type=$(echo $line | cut -d',' -f7)
replicate=$(echo $line | cut -d',' -f8)
rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
wd=$(echo $line | cut -d',' -f9 | cut -d':' -f1)
if [[ "$seq" = ATAC* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = Illumina ]]; then
      echo $srr "Align with STAR"
      star_atac
    elif [[ "$sequencer" = ABI ]]; then
      echo $srr "Align with BOWTIE"
      abi-align
    fi 
  done
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  investigator=$(echo $line | cut -d',' -f6)
  sequencer=$(echo $line | cut -d',' -f7)
  type=$(echo $line | cut -d',' -f8)
  replicate=$(echo $line | cut -d',' -f9)
  rep=$(echo $line | cut -d',' -f9 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f9 | cut -d'_' -f2)
  wd=$(echo $line | cut -d',' -f10 | cut -d':' -f1)
  cd $wd
  srr_line=$(echo $line | cut -d',' -f10 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = Illumina ]]; then
      echo $srr "Align with STAR"
      star_chip
    elif [[ "$sequencer" = ABI ]]; then
      echo $srr "Align with BOWTIE"
      abi-align
    fi 
  done
fi
if [[ "$seq" = RNA* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = Illumina ]]; then
      echo $srr "Align with STAR"
      star_rna
    elif [[ "$sequencer" = ABI ]]; then
      echo $srr "Align with BOWTIE"
      abi-align
    fi 
  done
fi
if [[ "$seq" = MNase* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = Illumina ]]; then
      echo $srr "Align with STAR"
      star_atac
    elif [[ "$sequencer" = ABI ]]; then
      echo $srr "Align with BOWTIE"
      abi-align
    fi 
  done
fi
done < $infodir/srr_files_${des}.txt
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
echo "#--------------------------------------------------------"
echo ""
echo ""; echo "Once files have been aligned, set duplicate = TRUE to mark and remove duplicates"; echo ""