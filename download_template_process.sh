#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory

##---------------------SETTING FUNCTIONS FOR PREPROCESSING--------------------------------
#Processes fastq files to remove adapters and low quality portions of reads
#Discards reads that are below 20bp
function adapter-trim {
  echo "TRIM_GALORE VERSION: `trim_galore --version`" >> $logdir/download_log_${des}.txt
  if [[ "$type" = "Paired" ]]; then
    echo "trim_galore -o `pwd` --${adapter} `pwd`/${srr}_1.fastq" >> $logdir/download_log_${des}.txt
    trim_galore --paired -o `pwd` --${adapter} `pwd`/${srr}_1.fastq `pwd`/${srr}_2.fastq
    mv ${srr}_1_trimmed.fq ${srr}_1.fastq
    mv ${srr}_2_trimmed.fq ${srr}_2.fastq
  elif [[ "$type" = "Single" ]]; then
    echo "trim_galore -o `pwd` --${adapter} `pwd`/${srr}.fastq" >> $logdir/download_log_${des}.txt
    trim_galore -o `pwd` --${adapter} `pwd`/${srr}.fastq
    mv ${srr}_trimmed.fq ${srr}.fastq
  fi
}

#Processes fastq files to mask low quality reads
function quality-trim {
  echo "Quality Trimming"
  echo "FASTX Toolkit (fastq_masker) version 0.0.14" >> $logdir/download_log_${des}.txt
  if [[ "$type" = "Paired" ]]; then
    echo "fastq_masker -q 20 -i `pwd`/${srr}_1.fastq -o `pwd`/${srr}_1_mask.fastq" >> $logdir/download_log_${des}.txt
    fastq_masker -q 20 -i `pwd`/${srr}_1.fastq -o `pwd`/${srr}_1_mask.fastq #Will mask low quality data
    mv ${srr}_1_mask.fastq ${srr}_1.fastq
    echo "fastq_masker -q 20 -i `pwd`/${srr}_2.fastq -o `pwd`/${srr}_2_mask.fastq" >> $logdir/download_log_${des}.txt
    fastq_masker -q 20 -i `pwd`/${srr}_2.fastq -o `pwd`/${srr}_mask_2.fastq #Will mask low quality data
    mv ${srr}_2_mask.fastq ${srr}_2.fastq
  elif [[ "$type" = "Single" ]]; then
    echo "fastq_masker -q 20 -i `pwd`/${srr}.fastq -o `pwd`/${srr}_mask.fastq" >> $logdir/download_log_${des}.txt
    fastq_masker -q 20 -i `pwd`/${srr}.fastq -o `pwd`/${srr}_mask.fastq #Will mask low quality data
    mv ${srr}_mask.fastq ${srr}.fastq
  fi
}

