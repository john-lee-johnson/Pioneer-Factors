#!/bin/bash

##---------------------SET ROTHENBERG BIOPROJECT NUMBER-----------------------------------
query=PRJNA146035

#Sets the working directories
infodir=/mnt/data1/John/Pioneer-Factors/info
maindir=/mnt/data1/John/Pioneer_Factors
scriptdir=/mnt/data1/John/Pioneer-Factors
filedir=/mnt/data1/John/Pioneer-Factors/sample_files
paralleldir=/mnt/data1/John/Pioneer-Factors/parallel_commands
r_dir=/mnt/data1/John/Pioneer-Factors/R_Scripts
datadir=/mnt/data1/VahediLab/PTF_Team/Data
dir0=$maindir

##---------------------SETTING FUNCTIONS--------------------------------------------------
#Creates parallel commands for fastq-dump
function fastq-download {
  fastq-dump -F -X 1 -Z ${srr} > $infodir/test.txt #Will download one line of the file to check for header information
    if [[ $(head $infodir/test.txt | sed -n 1p | cut -d":" -f1) = "@1" ]]; then
      echo "$srr" `pwd` >> $paralleldir/srr_download_rothenberg.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "$srr" `pwd` "-F" >> $paralleldir/srr_download_rothenberg.txt #If file does have proper header information, will keep the original header
    fi

echo `pwd` `pwd`"/${srr}.fastq illumina" >> $paralleldir/trim_galore_rothenberg.txt #Will remove sequencing adapters
echo "fastq_masker -q 20 -i `pwd`/${srr}_trimmed.fq -o `pwd`/${srr}_mask.fastq" >> $paralleldir/quality_masking_rothenberg.txt #Will mask low quality data
}

rm -f $infodir/srr_files_rothenberg.txt
rm -f $paralleldir/srr_download_rothenberg.txt
rm -f $paralleldir/trim_galore_rothenberg.txt
rm -f $paralleldir/quality_masking_rothenberg.txt

##---------------------DOWNLOAD BIOPROJECT SAMPLE DATA------------------------------------
esearch -db sra -query $query | efetch --format docsum | xtract -pattern DocumentSummary -element Title > $infodir/sample_files_rothenberg.txt
#Clean up sample names
sed -i -e '/FLDN/d' $infodir/sample_files_rothenberg.txt
sed -i -e 's/_/ /g' $infodir/sample_files_rothenberg.txt
sed -i -e 's/RNA-seq/RNA_seq/g' $infodir/sample_files_rothenberg.txt
sed -i -e '/ThyDN3 Input/d' $infodir/sample_files_rothenberg.txt
sed '' $infodir/sample_files_rothenberg.txt

##---------------------DOWNLOAD SRR NUMBERS DATA------------------------------------------
#Read in samples and generate a file of SRR numbers for Rothenberg Data
lines=$(wc -l $infodir/sample_files_rothenberg.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_rothenberg.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
seq=$(echo $line | cut -d':' -f2 | xargs | cut -d' ' -f2 | xargs)
cell=$(echo $line | cut -d':' -f2 | xargs | cut -d' ' -f1 | xargs)
rep=$(echo $line | cut -d':' -f2 | xargs | cut -d' ' -f3 | xargs)
echo $seq
if [[ "$seq" = RNA_seq* ]]; then
  cd $dir0/Data/RNA_seq ; mkdir -p $cell
  echo "$cell $seq" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_rothenberg.txt
else
  mark="$seq"
  if [[ "$seq" = GATA3* ]] || [[ "$seq" = PU.1 ]] || [[ "$seq" = Input ]]; then
    cd $dir0/Data/ChIP_seq/Transcription_Factor ; mkdir -p $cell ; cd $cell ; mkdir -p $mark ; cd $mark
    echo "$cell" "$mark" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_rothenberg.txt
  else
    cd $dir0/Data/ChIP_seq ; mkdir -p $cell ; cd $cell; mkdir -p $mark ; cd $mark
  echo "$cell $mark" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_rothenberg.txt
  fi
fi
done

##---------------------SRR DOWNLOAD DATA PARALLEL OUTPUT--------------------------------------
#Read in a file and set up parallel commands for later SRR download
while read line; do
seq=$(echo $line | cut -d' ' -f2)
cell=$(echo $line | cut -d' ' -f1)
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq/$cell
  srr_line=${line##"$cell $seq "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
  echo "$srr_line"
    echo "$srr $cell $seq"
    fastq-download
  done
else
  mark="$seq"
  if [[ "$seq" = GATA3* ]] || [[ "$seq" = PU.1 ]] || [[ "$seq" = Input ]]; then
    cd $dir0/Data/ChIP_seq/Transcription_Factor/$cell/$mark
    srr_line=${line##"$cell $seq "}
    IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo "$srr $cell $seq"
    fastq-download
  done
  else
    cd $dir0/Data/ChIP_seq/$cell/$mark
    srr_line=${line##"$cell $seq "}
    IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo "$srr $cell $seq"
    fastq-download
  done
  fi
fi
done < $infodir/srr_files_rothenberg.txt

