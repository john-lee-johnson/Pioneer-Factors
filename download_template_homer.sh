#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory
  
echo "MAKING HOMER TAG DIRECTORIES"; echo ""
if [ -f $paralleldir/homer_key_file.txt ]; then
  if [[ "$genome" = hg19 ]]; then
    cd $dir0/Analysis/Homer/Tag_Directories/; mkdir -p $des ; cd $des
    sed '' $paralleldir/homer_key_file.txt
    batchMakeTagDirectory.pl $paralleldir/homer_key_file.txt -cpu 40 -genome hg19 -format sam
  elif [[ "$genome" = mm10 ]]; then
    cd $dir0/Analysis/Homer/Tag_Directories/; mkdir -p $des ; cd $des
    sed '' $paralleldir/homer_key_file.txt
    batchMakeTagDirectory.pl $paralleldir/homer_key_file.txt -cpu 40 -genome mm10 -format sam    
  fi
else
  echo "Homer key file not found"
fi