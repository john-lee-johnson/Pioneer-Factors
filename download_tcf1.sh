#!/bin/bash
#This script will download sample information for Tcf1 and generate parallel commands
##---------------------SET BIOPROJECT NUMBER----------------------------------------------
query=PRJNA201421

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
      echo "$srr" `pwd` >> $paralleldir/srr_download_tcf1.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "$srr" `pwd` "-F" >> $paralleldir/srr_download_tcf1.txt #If file does have proper header information, will keep the original header
    fi
echo `pwd` `pwd`"/${srr}.fastq" >> $paralleldir/trim_galore_tcf1.txt #Will remove sequencing adapters
}

rm -f $infodir/srr_files_tcf1.txt
rm -f $paralleldir/srr_download_tcf1.txt
rm -f $paralleldir/trim_galore_tcf1.txt

##---------------------DOWNLOAD TCF1 BIOPROJECT SAMPLE DATA-------------------------------
esearch -db sra -query $query | efetch --format docsum | xtract -pattern DocumentSummary -element Title > $infodir/sample_files_tcf1.txt
#Clean up sample names
sed -i -nr '/Tc1f_WT|Input_Tcf1_WT/p' $infodir/sample_files_tcf1.txt
sed -i -e 's/Tc1f_WT/Tcf1_WT/g' $infodir/sample_files_tcf1.txt
sed -i -e 's/ChIP-Seq/ChIP_seq/g' $infodir/sample_files_tcf1.txt
sed -i -e 's/_WT/ Thy/g' $infodir/sample_files_tcf1.txt
sed -i -e 's/Input_Tcf1/Input/g' $infodir/sample_files_tcf1.txt
sed '' $infodir/sample_files_tcf1.txt

##---------------------DOWNLOAD SRR TCF1 DATA---------------------------------------------
#Read in samples and generate a file of SRR numbers for TCF1 Data
lines=$(wc -l $infodir/sample_files_tcf1.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_tcf1.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
seq=$(echo $line | cut -d':' -f2 | cut -d';' -f3 | xargs)
cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs)
echo $gsm
if [[ "$seq" = ChIP* ]]; then
  mark=$(echo $line | cut -d':' -f2 | xargs | cut -d';' -f1 | cut -d' ' -f1 | xargs)
  cell=$(echo $line | cut -d':' -f2 | xargs | cut -d';' -f1 | cut -d' ' -f2 | xargs)
  cd $dir0/Data/ChIP_seq/Transcription_Factor ; mkdir -p $cell ; cd $cell ; mkdir -p $mark ; cd $mark
  echo "$seq $mark $cell" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_tcf1.txt
fi
done


##---------------------SRR TCF1 DATA PARALLEL OUTPUT--------------------------------------
#Read in a file and set up parallel commands for later SRR download
while read line; do
seq=$(echo $line | cut -d' ' -f1)
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d' ' -f1)
  mark=$(echo $line | cut -d' ' -f2)
  cell=$(echo $line | cut -d' ' -f3)
  cd $dir0/Data/ChIP_seq/Transcription_Factor/$cell/$mark
  srr_line=${line##"$seq $mark $cell "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "download"
    fastq-download
  done
fi
done < $infodir/srr_files_tcf1.txt
