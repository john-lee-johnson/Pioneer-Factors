#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory

rm -f $infodir/sample_files_${des}_noreps.txt 
echo ""
echo "#----------Making New Sample File and Generating Homer Key File---------------------------"
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
  if [[ "$repnum" = 0 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
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
  if [[ "$repnum" = 0 ]]; then 
    filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    if echo "$mark" | grep -q Input; then
      echo "Input found, making Homer Tag Directory: `pwd`/${filename}.bam"
      echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
    else
      echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
    fi
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    if echo "$mark" | grep -q Input; then
      echo "Input found, making Homer Tag Directory"
      echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
    else
      echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
    fi
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  if [[ "$repnum" = 0 ]]; then 
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  echo "$line" >> $infodir/sample_files_${des}_noreps.txt
  echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  if [[ "$repnum" = 0 ]]; then 
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  echo "$line" >> $infodir/sample_files_${des}_noreps.txt
  echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
fi
done
echo ""
sed '' $infodir/sample_files_${des}_noreps.txt > $infodir/sample_files_${des}.txt 
