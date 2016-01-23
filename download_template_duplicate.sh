#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory

##---------------------SETTING FUNCTIONS FOR DUPLICATE REMOVAL----------------------------
#Duplicate removal should account for proper header information
#Generates commands for both marking and removal of duplicates

function mark_duplicates {
if [[ "$type" = Paired ]]; then
  if [[ $(head -1 ${srr}_1.fastq | cut -c 1-4) = "@SRR" ]]; then
    #If file does not have header information regarding location of the read on the lane, will not look for optical duplicates
    echo "PICARD VERSION: `java -jar $PICARD MarkDuplicates --version`" >> $logdir/download_log_${des}.txt
    echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_removed.bam READ_NAME_REGEX=null REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_${des}.txt
    echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_marked.bam READ_NAME_REGEX=null METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_${des}.txt
  else
    #If file has proper header information regarding location of the read on the lane, will look for optical duplicates
    echo "PICARD VERSION: `java -jar $PICARD MarkDuplicates --version`" >> $logdir/download_log_${des}.txt
    echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_removed.bam REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_${des}.txt
    echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_marked.bam METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_${des}.txt
  fi
elif [[ "$type" = Single ]]; then
  if [[ $(head -1 ${srr}.fastq | cut -c 1-4) = "@SRR" ]]; then
    #If file does not have header information regarding location of the read on the lane, will not look for optical duplicates
    echo "PICARD VERSION: `java -jar $PICARD MarkDuplicates --version`" >> $logdir/download_log_${des}.txt
    echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_removed.bam READ_NAME_REGEX=null REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_${des}.txt
    echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_marked.bam READ_NAME_REGEX=null METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_${des}.txt
  else
    #If file has proper header information regarding location of the read on the lane, will look for optical duplicates
    echo "PICARD VERSION: `java -jar $PICARD MarkDuplicates --version`" >> $logdir/download_log_${des}.txt
    echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_removed.bam REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_${des}.txt
    echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_marked.bam METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_${des}.txt
  fi
fi

}

function mark_duplicates_abi {
  #Will not look for optical duplicates because ABI SOLiD 
  echo "PICARD VERSION: `java -jar $PICARD MarkDuplicates --version`" >> $logdir/download_log_${des}.txt
  echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_removed.bam READ_NAME_REGEX=null REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_${des}.txt
  echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_marked.bam READ_NAME_REGEX=null METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_${des}.txt
}


echo "#--------------REMOVE DUPLICATES--------------------"
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
    echo $srr "duplicate removal"
    mkdir -p tmp
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates_abi
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates_abi
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
    echo $srr "duplicate removal"
    mkdir -p tmp
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates_abi
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates_abi
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
    echo $srr "duplicate removal"
    mkdir -p tmp
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates_abi
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates_abi
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
    echo $srr "duplicate removal"
    mkdir -p tmp
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates_abi
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates_abi
      fi
    fi
  done
fi
done < $infodir/srr_files_${des}.txt

if [[ $removed = true ]]; then
  if [ -f $paralleldir/remove_duplicates_${des}.txt ]; then
  echo "#--------------REMOVING DUPLICATES------------"
  echo "#--------------REMOVING DUPLICATES-------------------------" >> $logdir/download_log_${des}.txt
  parallel --xapply --dryrun -j 30 -- < $paralleldir/remove_duplicates_${des}.txt >> $logdir/download_log_${des}.txt
  echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
  parallel --xapply --dryrun -j 30 -- < $paralleldir/remove_duplicates_${des}.txt
  parallel --xapply -j 30 -- < $paralleldir/remove_duplicates_${des}.txt
fi
else
if [ -f $paralleldir/mark_duplicates_${des}.txt ]; then
  echo "Marking, but not removing duplicates"
  echo "#--------------MARKING DUPLICATES-----------------------------" >> $logdir/download_log_${des}.txt
  parallel --xapply --dryrun -j 30 -- < $paralleldir/mark_duplicates_${des}.txt >> $logdir/download_log_${des}.txt
  echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
  parallel --xapply --dryrun -j 30 -- < $paralleldir/mark_duplicates_${des}.txt
  parallel --xapply -j 30 -- < $paralleldir/mark_duplicates_${des}.txt
fi
fi

echo""; echo "#--------------COPY DEDUPLICATED FILES--------------------"
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
    if [[ $removed = true ]]; then
      echo $srr "Removing Duplicated Sequences"
      cp ${srr}_removed.bam ${srr}.bam
    else
      echo $srr "Keeping Duplicated Sequences"
      cp ${srr}_marked.bam ${srr}.bam
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
    if [[ $removed = true ]]; then
      echo $srr "Removing Duplicated Sequences"
      cp ${srr}_removed.bam ${srr}.bam
    else
      echo $srr "Keeping Duplicated Sequences"
      cp ${srr}_marked.bam ${srr}.bam
    fi
  done
fi
if [[ "$seq" = RNA* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ $removed = true ]]; then
      echo $srr "Removing Duplicated Sequences"
      cp ${srr}_removed.bam ${srr}.bam
    else
      echo $srr "Keeping Duplicated Sequences"
      cp ${srr}_marked.bam ${srr}.bam
    fi
  done
fi
if [[ "$seq" = MNase* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ $removed = true ]]; then
      echo $srr "Removing Duplicated Sequences"
      cp ${srr}_removed.bam ${srr}.bam
    else
      echo $srr "Keeping Duplicated Sequences"
      cp ${srr}_marked.bam ${srr}.bam
    fi
  done
fi
done < $infodir/srr_files_${des}.txt
echo "Now that duplicated sequences have been dealt with, set test4=true"