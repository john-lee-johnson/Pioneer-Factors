#!/bin/bash

##---------------------SET BIOPROJECT NUMBER----------------------------------------------
query=PRJNA218857

##---------------------SET SAMPLE DESCRIPTION---------------------------------------------
#Set des as a string that will identify the files being output in this script
des=mnase_bmdm

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
      echo "fastq-dump ${srr} -O ${wd}" >> $paralleldir/srr_download_${des}.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "fastq-dump -F ${srr} -O ${wd}"  >> $paralleldir/srr_download_${des}.txt #If file does have proper header information, will keep the original header
    fi
}

function fastq-download-paired {
  fastq-dump -F -X 1 -Z ${srr} > $infodir/test.txt #Will download one line of the file to check for header information
    echo "$(head $infodir/test.txt | sed -n 1p | cut -d":" -f1)"
    if [[ $(head $infodir/test.txt | sed -n 1p | cut -d":" -f1) = "@1" ]]; then
      echo "fastq-dump -I --split-files ${srr} -O ${wd}" >> $paralleldir/srr_download_${des}.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "fastq-dump -I --split-files -F ${srr} -O ${wd}"  >> $paralleldir/srr_download_${des}.txt #If file does have proper header information, will keep the original header
    fi
}

rm -f $infodir/srr_files_${des}.txt
rm -f $paralleldir/srr_download_${des}.txt
rm -f $paralleldir/trim_galore_${des}.txt
rm -f $paralleldir/homer_key_file.txt

##---------------------DOWNLOAD BIOPROJECT SAMPLE DATA------------------------------------
esearch -db sra -query $query | efetch --format docsum | xtract -pattern DocumentSummary -element Title > $infodir/sample_files_${des}.txt


##---------------------CLEAN UP SAMPLE FILE-----------------------------------------------
#The format for the sample file needs to be in the following format:
#GSM: Cell; Species; Seq; (Mark)
#Species should look like 'Mus musculus'
#Seq should look like 'ChIP_seq'
#If Seq is histone modification ChIP-seq, Seq should look like ChIP_seq; H3K4me1
#If Seq is transcription factor ChIP-seq, Seq should look like ChIP_seq; TCF1
#If Seq is input DNA ChIP-seq, Seq should look like ChIP_seq; Input

#Here are sample sed commands:
#Replace text with NewText
###sed -i -e 's/text/NewText/g' $infodir/sample_files_${des}.txt
#Delete lines with pattern
###sed -i -e '/pattern/d' $infodir/sample_files_${des}.txt
#Print lines matching pattern
###sed -i -n '/BMDM;/p' $infodir/sample_files_${des}.txt
sed -i -e 's/BMDM_UT_rep1/BMDM/g' $infodir/sample_files_${des}.txt
sed -i -e 's/BMDM_UT_rep2/BMDM/g' $infodir/sample_files_${des}.txt
sed -i -e 's/BMDM_UT_rep3/BMDM/g' $infodir/sample_files_${des}.txt
sed -i -e 's/BMDM_UT_rep4/BMDM/g' $infodir/sample_files_${des}.txt
sed -i -e 's/MNase-Seq/MNase_seq/g' $infodir/sample_files_${des}.txt
sed -i -n '/BMDM;/p' $infodir/sample_files_${des}.txt


echo "#--------------------SAMPLE FILES-----------------------"
sed '' $infodir/sample_files_${des}.txt
echo "#--------------------------------------------------------"

##---------------------DOWNLOAD SRR NUMBERS ----------------------------------------------
#Will read in GSM numbers from earlier and download the SRR numbers
#SRR File should have the following format:
#Example_seq Cell Directory SRR1 SRR2 SRR3 ... 
#Example_seq_Mark Cell Directory SRR1 SRR2 SRR3 ... 
echo ""
echo "#--------------------SRR FILES---------------------------"
lines=$(wc -l $infodir/sample_files_${des}.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs)
species=$(echo $line | cut -d':' -f2 | cut -d';' -f2 | xargs)
seq=$(echo $line | cut -d':' -f2 | cut -d';' -f3 | xargs)

if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; mkdir -p $cell ; cd $cell
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq to `pwd`"
  echo "$seq $cell `pwd`" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
  #esearch -db sra -query GSM1228683 | efetch --format runinfo | cut -d',' -f 16 | grep PAIRED
fi
if [[ "$seq" = ChIP* ]]; then
  mark=$(echo $line | cut -d':' -f2 | cut -d';' -f4 | xargs)
  cd $dir0/Data/ChIP_seq ; mkdir -p $cell/$mark ; cd $cell
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq $mark to `pwd`"
  echo "${seq}_${mark} $cell `pwd`" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; mkdir -p $cell ; cd $cell
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq to `pwd`"
  echo "$seq $cell `pwd`" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; mkdir -p $cell ; cd $cell
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq to `pwd`"
  echo "$seq $cell `pwd`" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
done
echo "#--------------------------------------------------------"
echo ""
echo "--------------------CHECK FILE FORMAT--------------------"
echo "SRR File should have the following format Example_seq_(Mark) Cell Directory SRR1 SRR2 SRR3 ..." 
sed '' $infodir/srr_files_${des}.txt
echo "#--------------------------------------------------------"
echo ""

###---------------------SRR DATA PARALLEL COMMANDS OUTPUT---------------------------------
#Read in a file and set up parallel commands for later SRR download
echo "--------------------FASTQ-DUMP--------------------"
while read line; do
seq=$(echo $line | cut -d' ' -f1)
cell=$(echo $line | cut -d' ' -f2)
wd=$(echo $line | cut -d' ' -f3)
if [[ "$seq" = ATAC* ]]; then
  cd $wd
  srr_line=${line##"$seq $cell $wd "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "download"
    fastq-download
  done
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d' ' -f1 | cut -d'_' -f1)
  mark=$(echo $line | cut -d' ' -f1 | cut -d'_' -f2)
  cd $wd
  srr_line=${line##"$seq_$mark $cell $wd "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "download"
    fastq-download
  done
fi
if [[ "$seq" = RNA* ]]
then
  cd $wd
  srr_line=${line##"$seq $cell $wd "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "download"
    fastq-download
  done
fi
if [[ "$seq" = MNase* ]]
then
  cd $wd
  srr_line=${line##"$seq $cell $wd "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "download"
    #fastq-download
    fastq-download-paired
  done
fi
done < $infodir/srr_files_${des}.txt
echo "#--------------------------------------------------------"
echo ""
echo "-----------------CHECK PARALLEL COMMANDS-----------------"
sed '' $paralleldir/srr_download_${des}.txt
echo "#--------------------------------------------------------"
echo ""

echo "-----------------SRR DOWNLOAD PARALLEL----------------"
parallel --dryrun --xapply -j 40 -- < $paralleldir/srr_download_${des}.txt
#parallel --xapply -j 40 -- < $paralleldir/srr_download_${des}.txt
#rm /mnt/data0/ncbi/public/sra/*
#fastq-dump -I --split-files -F SRR976078 -O /mnt/data1/John/Pioneer_Factors/Data/MNase_seq/BMDM
#rm /mnt/data0/ncbi/public/sra/*
#fastq-dump -I --split-files -F SRR976077 -O /mnt/data1/John/Pioneer_Factors/Data/MNase_seq/BMDM
#rm /mnt/data0/ncbi/public/sra/*
echo "#--------------------------------------------------------"
echo ""

###--------------COMBINE DUPLICATES--------------------------------
#If samples should be combined as duplicates that should be set here
#
echo "--------------------MNase Alignment--------------------"
echo "The following will be combined as duplicates"
cd $dir0/Data/MNase_seq/BMDM
#cp /mnt/data0/John/bowtie_mm10/* .
#bowtie mm10 -p 35 -v 3 -m 1 -S -I 0 -X 250 -1 SRR976077_1.fastq -2 SRR976077_2.fastq SRR976077.sam
#samtools view -h -b -o SRR976077.bam SRR976077.sam
#bowtie mm10 -p 35 -v 3 -m 1 -S -I 0 -X 250 -1 SRR976078_1.fastq -2 SRR976078_2.fastq SRR976078.sam
#samtools view -h -b -o SRR976078.bam SRR976078.sam
#bowtie mm10 -p 35 -v 3 -m 1 -S -I 0 -X 250 -1 SRR976079_1.fastq -2 SRR976079_2.fastq SRR976079.sam
#samtools view -h -b -o SRR976079.bam SRR976079.sam
#bowtie mm10 -p 35 -v 3 -m 1 -S -I 0 -X 250 -1 SRR976080_1.fastq -2 SRR976080_2.fastq SRR976080.sam
#samtools view -h -b -o SRR976080.bam SRR976080.sam
#bowtie -v 3 -m 1 -S -I 0 -X 250

#pigz -f -p 40 SRR976077_1.fastq
#pigz -f -p 40 SRR976077_2.fastq
#pigz -f -p 40 SRR976078_1.fastq
#pigz -f -p 40 SRR976078_2.fastq
#pigz -f -p 40 SRR976079_1.fastq
#pigz -f -p 40 SRR976079_2.fastq
#pigz -f -p 40 SRR976080_1.fastq
#pigz -f -p 40 SRR976080_2.fastq

samtools view -b -F 4 SRR976077.bam > SRR976077_mapped.bam
samtools view -b -F 4 SRR976078.bam > SRR976078_mapped.bam
samtools view -b -F 4 SRR976079.bam > SRR976079_mapped.bam
samtools view -b -F 4 SRR976080.bam > SRR976080_mapped.bam
