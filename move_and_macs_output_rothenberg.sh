#!/bin/bash

rm -f $paralleldir/macs_parallel_atac_srr*.txt
rm -f $paralleldir/macs_parallel_chip_srr*.txt
rm -f $paralleldir/homer_key_file.txt
rm -f $paralleldir/macs_parallel_atac.txt
rm -f $paralleldir/macs_parallel_chip.txt

#---------------------Setting Directories-------------------------------------------------
maindir=/mnt/data1/John/Pioneer_Factors
scriptdir=/mnt/data1/John/Pioneer-Factors
filedir=/mnt/data1/John/Pioneer-Factors/sample_files
paralleldir=/mnt/data1/John/Pioneer-Factors/parallel_commands
r_dir=/mnt/data1/John/Pioneer-Factors/R_Scripts
infodir=/mnt/data1/John/Pioneer-Factors/info
datadir=/mnt/data1/VahediLab/PTF_Team/Data
dir0=$maindir
cd ${dir0}

##---------------------COMBINE BAM FILES ROTHENBERG---------------------------------------
#Combine fastq files, compress to .gz, and align
awk '!a[$2,$3]++' $infodir/sample_files_rothenberg.txt > $infodir/sample_files_rothenberg_noreps.txt

lines=$(wc -l $infodir/sample_files_rothenberg_noreps.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_rothenberg_noreps.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
seq=$(echo $line | cut -d':' -f2 | xargs | cut -d' ' -f2 | xargs)
cell=$(echo $line | cut -d':' -f2 | xargs | cut -d' ' -f1 | xargs)
rep=$(echo $line | cut -d':' -f2 | xargs | cut -d' ' -f3 | xargs)
if [[ "$seq" = RNA_seq* ]]; then
  cd $dir0/Data/RNA_seq/$cell
  filename=${cell}_${seq}_Rothenberg_mm10_${gsm}
  samtools index ${filename}.bam
  cp ${filename}.bam $datadir/RNA_seq/bam
  cp ${filename}.bam.bai $datadir/RNA_seq/bam
  cp ${filename}.fastq.gz $datadir/RNA_seq/fastq
else
  mark="$seq"
  if [[ "$seq" = GATA3* ]] || [[ "$seq" = PU.1 ]] || [[ "$seq" = Input ]]; then
    cd $dir0/Data/ChIP_seq/Transcription_Factor/$cell/$mark
    filename=${cell}_ChIP_seq_${seq}_Rothenberg_mm10_${gsm}
    samtools index ${filename}.bam
    cp ${filename}.bam $datadir/ChIP_seq/bam
    cp ${filename}.bam.bai $datadir/ChIP_seq/bam
    cp ${filename}.fastq.gz $datadir/ChIP_seq/fastq
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  else
    cd $dir0/Data/ChIP_seq/$cell/$mark
    filename=${cell}_ChIP_seq_${seq}_Rothenberg_mm10_${gsm}
    samtools index ${filename}.bam
    cp ${filename}.bam $datadir/ChIP_seq/bam
    cp ${filename}.bam.bai $datadir/ChIP_seq/bam
    cp ${filename}.fastq.gz $datadir/ChIP_seq/fastq
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
fi
done

cd $dir0/Analysis/Homer/Tag_Directories/Rothenberg
batchMakeTagDirectory.pl $paralleldir/homer_key_file.txt -cpu 40 -genome mm10 -format sam
