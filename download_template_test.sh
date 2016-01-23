#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject

##---------------------SET BIOPROJECT NUMBER----------------------------------------------
query=SET_QUERY_NUMBER_HERE

##---------------------SET SAMPLE DESCRIPTION---------------------------------------------
#Set description string that will identify the info and parallel command files being output in this script
#Examples: mnase_bmdm, amit, mnase_dp, etc
des=SET_DESCRIPTION_HERE

#Sets the working directories
infodir=/mnt/data1/John/Pioneer-Factors/info
maindir=/mnt/data1/John/Pioneer_Factors
scriptdir=/mnt/data1/John/Pioneer-Factors
filedir=/mnt/data1/John/Pioneer-Factors/sample_files
paralleldir=/mnt/data1/John/Pioneer-Factors/parallel_commands
r_dir=/mnt/data1/John/Pioneer-Factors/R_Scripts
datadir=/mnt/data1/VahediLab/PTF_Team/Data
dir0=$maindir

rm -f $infodir/srr_files_${des}.txt
rm -f $paralleldir/srr_download_${des}.txt

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
    if [[ $(head $infodir/test.txt | sed -n 1p | cut -d":" -f1) = "@1" ]]; then
      echo "fastq-dump ${srr} -O ${wd}" >> $paralleldir/srr_download_${des}.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "fastq-dump -F ${srr} -O ${wd}"  >> $paralleldir/srr_download_${des}.txt #If file does have proper header information, will keep the original header
    fi
}
function abi-download {
  fastq-dump -F -X 1 -Z ${srr} > $infodir/test.txt #Will download one line of the file to check for header information
    if [[ $(head $infodir/test.txt | sed -n 1p | cut -d":" -f1) = "@1" ]]; then
      echo "fastq-dump ${srr} -O ${wd}" >> $paralleldir/srr_download_${des}.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "fastq-dump -F ${srr} -O ${wd}"  >> $paralleldir/srr_download_${des}.txt #If file does have proper header information, will keep the original header
    fi
}

##---------------------DOWNLOAD BIOPROJECT SAMPLE DATA------------------------------------
esearch -db sra -query $query | efetch --format docsum | xtract -pattern DocumentSummary -element Title > $infodir/sample_files_${des}.txt


##---------------------CLEAN UP SAMPLE FILE-----------------------------------------------
#The format for the sample file needs to be in the following format:
#GSM,Cell,Species,Seq,(Mark),Illumina/ABI,Paired/Single,Rep1
#Species should look like 'Mus musculus'
#Seq should look like 'ChIP_seq'
#If Seq is histone modification ChIP-seq, Seq should look like ChIP_seq,H3K4me1
#If Seq is transcription factor ChIP-seq, Seq should look like ChIP_seq,TCF1
#If Seq is input DNA ChIP-seq, Seq should look like ChIP_seq,Input
#If there are biological replicates, Rep should look like Rep_1 or Rep_2
#Here are sample sed commands:
#Replace text with NewText
###sed -i -e 's/text/NewText/g' $infodir/sample_files_${des}.txt
#Delete lines with pattern
###sed -i -e '/pattern/d' $infodir/sample_files_${des}.txt
#Print lines matching pattern
###sed -i -n '/pattern;/p' $infodir/sample_files_${des}.txt
#sed -i -e 's/text/NewText/g' $infodir/sample_files_${des}.txt
#sed -i -e 's/text/NewText/g' $infodir/sample_files_${des}.txt
#sed -i -e 's/text/NewText/g' $infodir/sample_files_${des}.txt
#sed -i -n '/pattern;/p' $infodir/sample_files_${des}.txt

echo "#--------------------SAMPLE FILES-----------------------"
sed '' $infodir/sample_files_${des}.txt
echo "#--------------------------------------------------------"

##---------------------DOWNLOAD SRR NUMBERS ----------------------------------------------
#Will read in GSM numbers from earlier and download the SRR numbers
#SRR File should have the following format:
#Example_seq (Mark) Cell Dir SRR1 SRR2 SRR3 ... 
echo ""
echo "#--------------------SRR FILES---------------------------"
lines=$(wc -l $infodir/sample_files_${des}.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}.txt)
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; mkdir -p $cell ; if [[ "$rep" = Rep ]]; then mkdir -p $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq at `pwd`"
  echo "$seq $cell `pwd`" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  rep=$(echo $line | cut -d',' -f8)
  cd $dir0/Data/ChIP_seq ; mkdir -p $cell/$mark ; cd $cell/$mark ; if [[ "$rep" = Rep ]]; then mkdir -p $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq $mark at `pwd`"
  echo "$seq $mark $cell `pwd`" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; mkdir -p $cell ; cd $cell ; if [[ "$rep" = Rep ]]; then mkdir -p $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq at `pwd`"
  echo "$seq $cell `pwd`" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; mkdir -p $cell ; cd $cell ; if [[ "$rep" = Rep ]]; then mkdir -p $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq at `pwd`"
  echo "$seq $cell `pwd`" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
done
echo "#--------------------------------------------------------"
echo "SRR File should have the following format Example_seq (Mark) Cell Dir SRR1 SRR2 SRR3 ..." 
sed '' $infodir/srr_files_${des}.txt

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
  cell=$(echo $line | cut -d' ' -f3)
  wd=$(echo $line | cut -d' ' -f4)
  seq=$(echo $line | cut -d' ' -f1)
  mark=$(echo $line | cut -d' ' -f2)
  cd $wd
  srr_line=${line##"$seq $mark $cell $wd "}
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
    fastq-download
  done
fi
done < $infodir/srr_files_${des}.txt
echo "#--------------------------------------------------------"
echo ""
echo "-----------------CHECK PARALLEL COMMANDS-----------------"
sed '' $paralleldir/srr_download_${des}.txt
echo "#--------------------------------------------------------"
echo ""