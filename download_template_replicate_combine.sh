#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory


rm -f $infodir/sample_files_${des}_noreps.txt
rm -f $paralleldir/merge_bam_${des}.txt

##-----------SETTING FUNCTIONS FOR COMBINING BIOLOGICAL REPLICATES-------------------------
#Combines bam files and fastq files for biological replicates (Illumina)
function replicate-combine-illumina {
unset input
unset reads
for a in $(ls *GSM*.bam | grep -v STAR | grep -v adapter | grep -v marked | grep -v removed | grep -v 'sort' | grep -v reads | grep -v combine); do
  input=$(echo $input" INPUT="`pwd`/$a)
done
echo "PICARD VERSION MergeSamFiles: `java -jar $PICARD MergeSamFiles --version`" >> $logdir/download_log_${des}.txt
echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MergeSamFiles $input ""OUTPUT="`pwd`/${filename}_combine.bam >> $paralleldir/merge_bam_${des}.txt
if [[ "$type" = Paired ]]; then
  echo "Skipping combining of fastq files because Paired-end"
else
  allFastq=$(ls *GSM*.fastq.gz | grep -v STAR | grep -v adapter | grep -v marked | grep -v removed | grep -v 'sort' | grep -v reads)
  cat $allFastq > ${filename}_combine.fastq.gz
fi  
}
#Combines bam files for biological replicates (ABI), but skips ABI fasta and qual files
function replicate-combine-abi {
unset input
unset reads
for a in $(ls *GSM*.bam | grep -v STAR | grep -v adapter | grep -v marked | grep -v removed | grep -v 'sort' | grep -v reads | grep -v combine); do
  input=$(echo $input" INPUT="`pwd`/$a)
done
echo "PICARD VERSION MergeSamFiles: `java -jar $PICARD MergeSamFiles --version`" >> $logdir/download_log_${des}.txt
echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MergeSamFiles $input ""OUTPUT="`pwd`/${filename}_combine.bam >> $paralleldir/merge_bam_${des}.txt
}


echo "#----------Making New Sample For Combining Biological Replicates---------------------------"
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
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
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
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  if [[ "$repnum" = 0 ]]; then 
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  echo "$line" >> $infodir/sample_files_${des}_noreps.txt
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
  fi
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  if [[ "$repnum" = 0 ]]; then 
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  echo "$line" >> $infodir/sample_files_${des}_noreps.txt
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
  fi
fi
done
sed '' $infodir/sample_files_${des}_noreps.txt
echo "#------------------------------------------------------------------"
echo "To combine biological replicates, set replicate_combine=true"
echo ""

echo "#-------------COMBINING BIOLOGICAL REPLICATES---------------------------"
lines=$(wc -l $infodir/sample_files_${des}_noreps.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}_noreps.txt)
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
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    if [[ "$sequencer" = Illumina ]]; then
      replicate-combine-illumina
    elif [[ "$sequencer" = ABI ]]; then
      replicate-combine-abi
    fi 
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
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
    if [[ "$sequencer" = Illumina ]]; then
      replicate-combine-illumina
    elif [[ "$sequencer" = ABI ]]; then
      replicate-combine-abi
    fi 
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    if [[ "$sequencer" = Illumina ]]; then
      replicate-combine-illumina
    elif [[ "$sequencer" = ABI ]]; then
      replicate-combine-abi
    fi 
  fi
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    if [[ "$sequencer" = Illumina ]]; then
      replicate-combine-illumina
    elif [[ "$sequencer" = ABI ]]; then
      replicate-combine-abi
    fi 
  fi
fi
done
echo "#--------------------------------------------------------"
echo ""
sed '' $infodir/sample_files_${des}_noreps.txt

echo "COMBINING BAM FILES"; echo ""
if [ -f $paralleldir/merge_bam_${des}.txt ]; then
echo "#---------------COMBINING BIO REPS--------------------------" >> $logdir/download_log_${des}.txt
parallel --xapply --dryrun -j 30 -- < $paralleldir/merge_bam_${des}.txt >> $logdir/download_log_${des}.txt
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
parallel --xapply --dryrun -j 30 -- < $paralleldir/merge_bam_${des}.txt
parallel --xapply -j 30 -- < $paralleldir/merge_bam_${des}.txt
echo "#-------------Clean Up Directories---------------------------"
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
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  if [[ "$repnum" -ne 0 ]]; then
    if [[ "$repnum" -eq 1 ]]; then 
      if [[ "$sequencer" = Illumina ]]; then
        if [[ "$type" = Paired ]]; then
          rm -f ${filename}.bam
          mv ${filename}_combine.bam ${filename}.bam
        elif [[ "$type" = Single ]]; then
          rm -f ${filename}.bam
          rm -f ${filename}.fastq.gz
          mv ${filename}_combine.bam ${filename}.bam
          mv ${filename}_combine.fastq.gz ${filename}.fastq.gz
        fi 
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
        mv ${filename}_combine.bam ${filename}.bam
      fi 
    else
      if [[ "$sequencer" = Illumina ]]; then
        rm -f ${filename}.bam
        rm -f ${filename}.fastq.gz
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
      fi 
    fi  
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
  filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
  if [[ "$repnum" -ne 0 ]]; then
    if [[ "$repnum" -eq 1 ]]; then 
      if [[ "$sequencer" = Illumina ]]; then
        if [[ "$type" = Paired ]]; then
          rm -f ${filename}.bam
          mv ${filename}_combine.bam ${filename}.bam
        elif [[ "$type" = Single ]]; then
          rm -f ${filename}.bam
          rm -f ${filename}.fastq.gz
          mv ${filename}_combine.bam ${filename}.bam
          mv ${filename}_combine.fastq.gz ${filename}.fastq.gz
        fi 
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
        mv ${filename}_combine.bam ${filename}.bam
      fi 
    else
      if [[ "$sequencer" = Illumina ]]; then
        rm -f ${filename}.bam
        rm -f ${filename}.fastq.gz
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
      fi 
    fi  
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  if [[ "$repnum" -ne 0 ]]; then
    if [[ "$repnum" -eq 1 ]]; then 
      if [[ "$sequencer" = Illumina ]]; then
        if [[ "$type" = Paired ]]; then
          rm -f ${filename}.bam
          mv ${filename}_combine.bam ${filename}.bam
        elif [[ "$type" = Single ]]; then
          rm -f ${filename}.bam
          rm -f ${filename}.fastq.gz
          mv ${filename}_combine.bam ${filename}.bam
          mv ${filename}_combine.fastq.gz ${filename}.fastq.gz
        fi 
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
        mv ${filename}_combine.bam ${filename}.bam
      fi 
    else
      if [[ "$sequencer" = Illumina ]]; then
        rm -f ${filename}.bam
        rm -f ${filename}.fastq.gz
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
      fi 
    fi  
  fi
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  if [[ "$repnum" -ne 0 ]]; then
    if [[ "$repnum" -eq 1 ]]; then 
      if [[ "$sequencer" = Illumina ]]; then
        if [[ "$type" = Paired ]]; then
          rm -f ${filename}.bam
          mv ${filename}_combine.bam ${filename}.bam
        elif [[ "$type" = Single ]]; then
          rm -f ${filename}.bam
          rm -f ${filename}.fastq.gz
          mv ${filename}_combine.bam ${filename}.bam
          mv ${filename}_combine.fastq.gz ${filename}.fastq.gz
        fi 
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
        mv ${filename}_combine.bam ${filename}.bam
      fi 
    else
      if [[ "$sequencer" = Illumina ]]; then
        rm -f ${filename}.bam
        rm -f ${filename}.fastq.gz
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
      fi 
    fi  
  fi
fi
done
fi