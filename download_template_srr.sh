#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory

##---------------------DOWNLOAD SRR NUMBERS ----------------------------------------------
#Will read in GSM numbers from SAMPLE FILE and download associated SRR numbers
#SRR File should have the following format:
#gsm,cell,species,seq,(mark),investigator,sequencer,type,replicate,wd:SRR1 SRR2 SRR3 ...
echo ""
echo "#--------------DOWNLOAD SRR NUMBERS---------------------------"
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
  cd $dir0/Data/ATAC_seq ; mkdir -p $cell/$investigator ; cd $cell/$investigator ; if [[ "$repnum" -ne 0 ]]; then mkdir -p $replicate; cd $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq at `pwd`"
  echo "$gsm,$cell,$species,$seq,$investigator,$sequencer,$type,$replicate,`pwd`:"`esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  investigator=$(echo $line | cut -d',' -f6)
  sequencer=$(echo $line | cut -d',' -f7)
  type=$(echo $line | cut -d',' -f8)
  replicate=$(echo $line | cut -d',' -f9)
  rep=$(echo $line | cut -d',' -f9 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f9| cut -d'_' -f2)
  cd $dir0/Data/ChIP_seq ; mkdir -p $cell/$mark/$investigator ; cd $cell/$mark/$investigator ; if [[ "$repnum" -ne 0 ]]; then mkdir -p $replicate; cd $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq $mark at `pwd`"
  echo "$gsm,$cell,$species,$seq,$mark,$investigator,$sequencer,$type,$replicate,`pwd`:"`esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; mkdir -p $cell/$investigator ; cd $cell/$investigator ; if [[ "$repnum" -ne 0 ]]; then mkdir -p $replicate; cd $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq at `pwd`"
  echo "$gsm,$cell,$species,$seq,$investigator,$sequencer,$type,$replicate,`pwd`:"`esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; mkdir -p $cell/$investigator ; cd $cell/$investigator ; if [[ "$repnum" -ne 0 ]]; then mkdir -p $replicate; cd $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq at `pwd`"
  echo "$gsm,$cell,$species,$seq,$investigator,$sequencer,$type,$replicate,`pwd`:"`esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
done
echo "--------------------------------------------------------"
echo "SRR File should have the following format gsm,cell,species,seq,(mark),investigator,sequencer,type,pwd:SRR1 SRR2 SRR3 ..." 
echo "";echo "#--------------SRR FILES---------------------------------"
sed '' $infodir/srr_files_${des}.txt
echo ""
