#!/bin/bash

#!/bin/bash

##---------------------SET THE BIOPROJECT NUMBER-----------------------------------------
query=PRJNA175500

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
#if [ ! -f ${srr}.fastq ]; then
  fastq-dump -F -X 1 -Z ${srr} > $infodir/test.txt #Will download one line of the file to check for header information
    if [[ $(head $infodir/test.txt | sed -n 1p | cut -d":" -f1) = "@1" ]]; then
      echo "$srr" `pwd` >> $paralleldir/srr_download_mnase.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "$srr" `pwd` "-F" >> $paralleldir/srr_download_mnase.txt #If file does have proper header information, will keep the original header
    fi
#else
#  echo "Raw fastq file exists!"
#fi
if [[ "$seq" = ATAC* ]]; then
      echo `pwd` `pwd`"/${srr}.fastq nextera" >> $paralleldir/trim_galore.txt #Will remove sequencing adapters
fi
if [[ "$seq" = RNA* ]]; then
      echo `pwd` `pwd`"/${srr}.fastq illumina" >> $paralleldir/trim_galore.txt #Will remove sequencing adapters
      #echo "fastq_masker -q 20 -i `pwd`/${srr}_trimmed.fq -o `pwd`/${srr}_mask.fastq" >> $paralleldir/quality_masking_rothenberg.txt
fi
}

rm -f $infodir/srr_files_mnase.txt
rm -f $paralleldir/srr_download_mnase.txt
rm -f $paralleldir/trim_galore_mnase.txt

##---------------------DOWNLOAD AMIT BIOPROJECT SAMPLE DATA-------------------------------
esearch -db sra -query $query | efetch --format docsum | xtract -pattern DocumentSummary -element Title > $infodir/sample_files_mnase.txt
#Clean up sample names
sed -i -e '/GSM1004653/!d' $infodir/sample_files_mnase.txt
sed -i -e 's/ merged mononucleosomes//g' $infodir/sample_files_mnase.txt
sed '' $infodir/sample_files_mnase.txt
