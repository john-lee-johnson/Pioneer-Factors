#!/bin/bash
#This script will download the Sample information from NCBI for Ido Amit data
#Generates scripts used for parallel download and trimming of adapters

##---------------------SET AMIT BIOPROJECT NUMBER-----------------------------------------
query=PRJNA257488

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
      echo "$srr" `pwd` >> $paralleldir/srr_download.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "$srr" `pwd` "-F" >> $paralleldir/srr_download.txt #If file does have proper header information, will keep the original header
    fi
if [[ "$seq" = ATAC* ]]; then
      echo `pwd` `pwd`"/${srr}.fastq nextera" >> $paralleldir/trim_galore.txt #Will remove sequencing adapters for ATAC-seq files
fi
if [[ "$seq" = RNA* ]]; then
      echo `pwd` `pwd`"/${srr}.fastq illumina" >> $paralleldir/trim_galore.txt #Will remove sequencing adapters from RNA-seq files
fi
}

rm -f $infodir/srr_files_amit.txt
rm -f $paralleldir/srr_download.txt
rm -f $paralleldir/trim_galore.txt

##---------------------DOWNLOAD BIOPROJECT SAMPLE DATA------------------------------------
esearch -db sra -query PRJNA257488 | efetch --format docsum | xtract -pattern DocumentSummary -element Title > $infodir/sample_files_amit.txt
#Clean up sample names
sed -i -e 's/Granulocytes/GN/g' $infodir/sample_files_amit.txt
sed -i -e 's/Granulocyte/GN/g' $infodir/sample_files_amit.txt
sed -i -e 's/Monocytes/Mono/g' $infodir/sample_files_amit.txt
sed -i -e 's/LT-HSC/LTHSC/g' $infodir/sample_files_amit.txt
sed -i -e 's/ST-HSC/STHSC/g' $infodir/sample_files_amit.txt
sed -i -e 's/LT_HSC/LTHSC/g' $infodir/sample_files_amit.txt
sed -i -e 's/ST_HSC/STHSC/g' $infodir/sample_files_amit.txt
sed -i -e  's/\<HSC\>/STHSC/g' $infodir/sample_files_amit.txt
sed -i -e 's/.ucsc//g' $infodir/sample_files_amit.txt
sed -i -e 's/ ATAC_seq//g' $infodir/sample_files_amit.txt
sed -i -e 's/ChIP-Seq/ChIP_seq/g' $infodir/sample_files_amit.txt
sed -i -e 's/RNA-Seq/RNA_seq/g' $infodir/sample_files_amit.txt
sed -i -e 's/OTHER/ATAC_seq/g' $infodir/sample_files_amit.txt

##---------------------DOWNLOAD SRR NUMBERS-----------------------------------------------
#Read in samples and generate a file of SRR numbers
lines=$(wc -l $infodir/sample_files_amit.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_amit.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
seq=$(echo $line | cut -d':' -f2 | cut -d';' -f3 | xargs)
cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs)
echo $gsm
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; mkdir -p $cell ; cd $cell
  echo "$seq $cell" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_amit.txt
fi
if [[ "$seq" = ChIP* ]]; then
  mark=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | cut -d'_' -f1 | xargs)
  cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | cut -d'_' -f2 | xargs)
  cd $dir0/Data/ChIP_seq ; mkdir -p $cell ; cd $cell ; mkdir -p $mark ; cd $mark
  echo "$seq $mark $cell" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_amit.txt
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; mkdir -p $cell ; cd $cell
  echo "$seq $cell" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_amit.txt
fi
done

##---------------------SRR DATA DOWNLOAD PARALLEL OUTPUT----------------------------------
#Read in a file and set up parallel commands for later SRR download
while read line; do
seq=$(echo $line | cut -d' ' -f1)
if [[ "$seq" = ATAC* ]]; then
  cell=$(echo $line | cut -d' ' -f2)
  cd $dir0/Data/ATAC_seq/$cell
  srr_line=${line##"$seq $cell "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "download"
    fastq-download
  done
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d' ' -f1)
  mark=$(echo $line | cut -d' ' -f2)
  cell=$(echo $line | cut -d' ' -f3)
  cd $dir0/Data/ChIP_seq/$cell/$mark
  srr_line=${line##"$seq $mark $cell "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "download"
    fastq-download
  done
fi
if [[ "$seq" = RNA* ]]
then
  cell=$(echo $line | cut -d' ' -f2)
  cd $dir0/Data/RNA_seq/$cell
  srr_line=${line##"$seq $cell "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "download"
    fastq-download
  done
fi
done < $infodir/srr_files_amit.txt
