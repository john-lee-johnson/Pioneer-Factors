#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory

echo "Combining Biological Replicates"; echo ""
echo "#--------------------Combining Biological Replicates---------------------------"
echo "#Part 1 Copying Files"
lines=$(wc -l $infodir/sample_files_${des}.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}.txt)
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
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell/$investigator
  if [[ "$repnum" -ne 0 ]]; then 
    cd $replicate
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    cp ${filename}.bam ..
    if [[ "$sequencer" = Illumina ]]; then
      if [[ "$type" = Paired ]]; then
        echo "Skipping copying of FASTQ files because Paired-End Reads"
      elif [[ "$type" = Single ]]; then
        echo "Copying Rep $repnum fastq to upper directory"
        cp ${filename}.fastq.gz ..  
      fi      
    elif [[ "$sequencer" = ABI ]]; then
      echo "Skipping ${filename} because ABI"
    fi 
  fi
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
  cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator
  if [[ "$repnum" -ne 0 ]]; then
    cd $replicate
    filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
    cp ${filename}.bam ..
    if [[ "$sequencer" = Illumina ]]; then
      if [[ "$type" = Paired ]]; then
        echo "Skipping copying of FASTQ files because Paired-End Reads"
      elif [[ "$type" = Single ]]; then
        echo "Copying Rep $repnum fastq to upper directory"
        cp ${filename}.fastq.gz ..  
      fi 
    elif [[ "$sequencer" = ABI ]]; then
      echo "Skipping ${filename} because ABI"
    fi 
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  if [[ "$repnum" -ne 0 ]]; then
    cd $replicate
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    cp ${filename}.bam ..
    if [[ "$sequencer" = Illumina ]]; then
      if [[ "$type" = Paired ]]; then
        echo "Skipping copying of FASTQ files because Paired-End Reads"
      elif [[ "$type" = Single ]]; then
        echo "Copying Rep $repnum fastq to upper directory"
        cp ${filename}.fastq.gz ..  
      fi 
    elif [[ "$sequencer" = ABI ]]; then
      echo "Skipping ${filename} because ABI"
    fi 
  fi
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  if [[ "$repnum" -ne 0 ]]; then
    cd $replicate
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    cp ${filename}.bam ..
    if [[ "$sequencer" = Illumina ]]; then
      if [[ "$type" = Paired ]]; then
        echo "Skipping copying of FASTQ files because Paired-End Reads"
      elif [[ "$type" = Single ]]; then
        echo "Copying Rep $repnum fastq to upper directory"
        cp ${filename}.fastq.gz ..  
      fi 
    elif [[ "$sequencer" = ABI ]]; then
      echo "Skipping ${filename} because ABI"
    fi 
  fi
fi
done