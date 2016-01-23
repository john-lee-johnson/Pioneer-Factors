#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory

##---------------------SETTING FUNCTIONS FOR DOWNLOAD-------------------------------------
#Download functions should account for proper header information and paired vs. single end data

#Creates parallel commands for fastq-dump single-end Illumina data
function fastq-download {
  fastq-dump -F -X 1 -Z ${srr} > $infodir/test.txt #Will download one line of the file to check for header information
    if [[ $(head $infodir/test.txt | sed -n 1p | cut -d":" -f1) = "@1" ]]; then
      echo "fastq-dump ${srr} -O ${wd}" >> $paralleldir/srr_download_${des}.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "fastq-dump -F ${srr} -O ${wd}"  >> $paralleldir/srr_download_${des}.txt #If file does have proper header information, will keep the original header
    fi
}
#Creates parallel commands for fastq-dump paired-end Illumina data
function fastq-download-paired {
  fastq-dump -F -X 1 -Z ${srr} > $infodir/test.txt #Will download one line of the file to check for header information
    if [[ $(head $infodir/test.txt | sed -n 1p | cut -d":" -f1) = "@1" ]]; then
      echo "fastq-dump -I --split-files ${srr} -O ${wd}" >> $paralleldir/srr_download_${des}.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "fastq-dump -I --split-files -F ${srr} -O ${wd}"  >> $paralleldir/srr_download_${des}.txt #If file does have proper header information, will keep the original header
    fi
}
#Creates parallel commands for fastq-dump single-end ABI data
function abi-download {
    echo "abi-dump -F ${srr} -O ${wd}" >> $paralleldir/srr_download_${des}.txt #If file does not have proper header information, will keep the SRR header
}
#Creates parallel commands for fastq-dump paired-end ABI data
function abi-download-paired {
    echo "abi-dump -F ${srr} -O ${wd}" >> $paralleldir/srr_download_${des}.txt #If file does not have proper header information, will keep the SRR header
}

#Checks if the SRA file is paired or single-end 
function check_sra_paired {
  sra_paired() {
    local SRA="$srr"
    local x=$(
      fastq-dump -I -X 1 -Z --split-spot "$SRA" 2>/dev/null \
        | awk '{if(NR % 2 == 1) print substr($1,length($1),1)}' \
        | uniq \
        | wc -l
    )
    [[ $x == 2 ]]
  }
if sra_paired "$1"; then
  echo "Downloading Paired-end Illumina data for $srr"
  fastq-download-paired
else
  echo "Downloading Single-end Illumina data for $srr"
  fastq-download
fi
}

echo "#--------------SETTING UP COMMANDS FOR DOWNLOAD-------------------"

###---------------------SRA TOOLKIT PARALLEL COMMANDS OUTPUT------------------------------
#Read in the SRR FILE, and set up parallel commands for later SRR download
while read line; do
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
wd=$(echo $line | cut -d',' -f9 | cut -d':' -f1)
if [[ "$seq" = ATAC* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        check_sra_paired
      elif [[ "$type" = "Single" ]]; then
        check_sra_paired
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "Paired-end ABI data"
        abi-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "Download Single-end ABI data"
        abi-download
      fi
    fi
  done
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
  wd=$(echo $line | cut -d',' -f10 | cut -d':' -f1)
  cd $wd
  srr_line=$(echo $line | cut -d',' -f10 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        check_sra_paired
      elif [[ "$type" = "Single" ]]; then
        check_sra_paired
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "Paired-end ABI data"
        abi-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "Download Single-end ABI data"
        abi-download
      fi
    fi
  done
fi
if [[ "$seq" = RNA* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        check_sra_paired
      elif [[ "$type" = "Single" ]]; then
        check_sra_paired
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "Paired-end ABI data"
        abi-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "Download Single-end ABI data"
        abi-download
      fi
    fi
  done
fi
if [[ "$seq" = MNase* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        check_sra_paired
      elif [[ "$type" = "Single" ]]; then
        check_sra_paired
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "Paired-end ABI data"
        abi-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "Download Single-end ABI data"
        abi-download
      fi
    fi
  done
fi
done < $infodir/srr_files_${des}.txt
echo "#--------------------------------------------------------"
echo ""
echo "#--------------CHECK PARALLEL COMMANDS-----------------"
sed '' $paralleldir/srr_download_${des}.txt
echo "#--------------------------------------------------------"
echo "Check to make sure download commands are correct before setting download = true"
echo ""
