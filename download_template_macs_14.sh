#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory

##--------------MACS 14 OUTPUT COMMANDS COMBINED BAM FILE-------------------------------
#Read in a file and output commands to a file for MACS peak calling 14
lines=$(wc -l $infodir/sample_files_${des}.txt  | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}.txt )
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
  cd $dir0/Data/ATAC_seq ; cd $cell/$investigator ; mkdir -p macs
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  if [[ "$type" = Paired ]]; then
    fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/$des/$filename/tagInfo.txt" | cut -d"=" -f2 | xargs)
    echo "`pwd`/${filename}.bam `pwd`/macs/${filename}" "$fragLength">> $paralleldir/macs_parallel_${des}.txt  
  elif [[ "$type" = Single ]]; then
    echo "`pwd`/${filename}.bam `pwd`/macs/${filename}">> $paralleldir/macs_parallel_${des}.txt  
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
  filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
  if echo "$mark" | grep -q Input; then
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator ; mkdir -p macs
    echo "Input found: `pwd`/${filename}.bam"
    if [[ "$type" = Paired ]]; then
      fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/$des/$filename/tagInfo.txt" | cut -d"=" -f2 | xargs)
      echo "`pwd`/${filename}.bam `pwd`/macs/${filename}" "$fragLength">> $paralleldir/macs_parallel_${des}.txt  
    elif [[ "$type" = Single ]]; then
      echo "`pwd`/${filename}.bam `pwd`/macs/${filename}">> $paralleldir/macs_parallel_${des}.txt  
    fi
  else
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator ; mkdir -p macs
    if [[ "$type" = Paired ]]; then
      fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/$des/$filename/tagInfo.txt" | cut -d"=" -f2 | xargs)
      echo "`pwd`/${filename}.bam `pwd`/macs/${filename}" "$fragLength">> $paralleldir/macs_parallel_${des}.txt  
    elif [[ "$type" = Single ]]; then
      echo "`pwd`/${filename}.bam `pwd`/macs/${filename}">> $paralleldir/macs_parallel_${des}.txt  
    fi
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
done