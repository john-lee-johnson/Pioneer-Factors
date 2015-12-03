#!/bin/bash
##This script will read the SRR files from srr_files.txt, download and convert them to fastq, align them, output to BAM files, and set up commands for MACS peak calling
##This is the first script in the series to be run

#Sets the working directories
dir0=/mnt/data1/John/Pioneer_Factors
filedir=/mnt/data1/John/Pioneer-Factors/sample_files
datadir=/mnt/data1/VahediLab/PTF_Team/Data
paralleldir=/mnt/data1/John/Pioneer-Factors/parallel_commands

##---------------------SETTING FUNCTIONS--------------------------------------------------
#STAR generates genome index
function star_generate {
STAR --runMode genomeGenerate --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --genomeFastaFiles /mnt/data0/John/GRCm38p4_mm10_pa_only.fa --sjdbGTFfile /mnt/data0/John/gencode.vM6.annotation.gtf
}

#STAR alignment for ATAC-seq
function star_atac {
STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate 
mv Aligned.sortedByCoord.out.bam STAR_${srr}.bam #Keeps the bam file that STAR generates
mv Log.final.out ${srr}_Log.final.out #Keeps a copy of the aligner log
STAR --runMode inputAlignmentsFromBAM --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --inputBAMfile STAR_${srr}.bam --bamRemoveDuplicatesType UniqueIdentical
mv Processed.out.bam dedupSTAR_${srr}.bam #Keeps a bam file with PCR duplicates removed via STAR
}

#STAR alignment for ChIP-seq
function star_chip {
STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate 
mv Aligned.sortedByCoord.out.bam STAR_${srr}.bam #Keeps the bam file that STAR generates
mv Log.final.out ${srr}_Log.final.out #Keeps a copy of the aligner log
STAR --runMode inputAlignmentsFromBAM --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --inputBAMfile STAR_${srr}.bam --bamRemoveDuplicatesType UniqueIdentical
mv Processed.out.bam dedup_${srr}.bam #Keeps a bam file with PCR duplicates removed via STAR
}

function mark_duplicates {
if [[ $(head -1 ${srr}.fastq | cut -c 1-4) = "@SRR" ]]; then
  #If file does not have header information regarding location of the read on the lane, will not look for optical duplicates
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/STAR_${srr}.bam OUTPUT="`pwd`"/${srr}.bam READ_NAME_REGEX=null REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_rothenberg.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/STAR_${srr}.bam OUTPUT="`pwd`"/dupsMarked_${srr}.bam READ_NAME_REGEX=null METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_rothenberg.txt

else
  #If file has proper header information regarding location of the read on the lane, will look for optical duplicates
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/STAR_${srr}.bam OUTPUT="`pwd`"/${srr}.bam REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_rothenberg.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/STAR_${srr}.bam OUTPUT="`pwd`"/dupsMarked_${srr}.bam METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_rothenberg.txt
fi
}

#STAR alignment for RNA-seq 2-pass
function star_rna {
STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15
mv Aligned.sortedByCoord.out.bam ${srr}.bam #Keeps the bam file that STAR generates
mv Log.final.out ${srr}_Log.final.out #Keeps a copy of the aligner log
STAR --runMode inputAlignmentsFromBAM --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --inputBAMfile ${srr}.bam --bamRemoveDuplicatesType UniqueIdentical
mv Processed.out.bam dedupSTAR_${srr}.bam #Keeps a bam file with PCR duplicates removed via STAR
}

#STAR remove genome from memory
function star_remove {
STAR --genomeLoad Remove --genomeDir /mnt/data0/John/genome_GRCm38p4_M6
}

#Create correlation heatmaps between replicates
function heatmap {
for a in $(ls SRR*.bam | grep -v STAR); do
  samtools index $a
done
if [[ $(ls SRR*.bam | grep -v STAR | wc -l) > 1 ]]; then
unset input
for b in $(ls SRR*.bam | grep -v STAR); do
  input=$(echo $input" "`pwd`/$b)
done
  echo "bamCorrelate bins --bamfiles $input --corMethod spearman -o" `pwd`/"spearmanHeatmap.png" >> $paralleldir/heat_map_rothenberg.txt
  echo "bamCorrelate bins --bamfiles $input --corMethod pearson -o" `pwd`"/pearsonHeatmap.png" >> $paralleldir/heat_map_rothenberg.txt
fi
}

#Combine replicates
function bam_combine {
unset input
unset reads
for a in $(ls SRR*.bam | grep -v STAR); do
  input=$(echo $input" INPUT="`pwd`/$a)
  #reads=reads+$(samtools view -c $a)
done
echo "java -Xmx2g -jar $PICARD MergeSamFiles $input ""OUTPUT="`pwd`/${filename}.bam >> $paralleldir/merge_bam_rothenberg.txt
#echo $seq $cell "Dedup" $(samtools view -c ${filename}.bam) "Total" $reads >> $infodir/pcr_duplicates.txt
}

#Make sure genome is unloaded on first-run
cd $dir0
star_remove
#rm -f $infodir/pcr_duplicates.txt
rm -f $paralleldir/mark_duplicates_rothenberg.txt
rm -f $paralleldir/remove_duplicates_rothenberg.txt
rm -f $paralleldir/merge_bam_rothenberg.txt
rm -f $paralleldir/heat_map_rothenberg.txt

##---------------------DOWNLOAD SRR ROTHENBERG DATA---------------------------------------
#Download the SRR files for Amit Data
if [ -f $paralleldir/srr_download_rothenberg.txt ]; then
parallel --xapply --dryrun -j 40 --colsep ' ' -a $paralleldir/srr_download_rothenberg.txt  "fastq-dump {3} {1} -O {2}"
#parallel --xapply -j 40 --colsep ' ' -a $paralleldir/srr_download_rothenberg.txt  "fastq-dump {3} {1} -O {2}"
fi

##---------------------REMOVE ADAPTERS ROTHENBERG DATA------------------------------------
#Remove sequencing adapters from Amit Data
if [ -f $paralleldir/trim_galore_rothenberg.txt ]; then
parallel --xapply --dryrun -j 40 --colsep ' ' -a $paralleldir/trim_galore_rothenberg.txt "trim_galore -o {1} --{3} {2} --length 15 -q 20" #Removes Nextera Transposase Adapters
#parallel --xapply -j 40 --colsep ' ' -a $paralleldir/trim_galore_rothenberg.txt "trim_galore -o {1} --{3} {2} --length 15 -q 20"
fi
#parallel --xapply --dryrun -j 40 -- < $paralleldir/quality_masking_rothenberg.txt
#parallel --xapply -j 40 -- < $paralleldir/quality_masking_rothenberg.txt

##---------------------SRR ROTHENBERG DATA PARALLEL OUTPUT--------------------------------
#Read in a file and set up parallel commands for later SRR download
while read line; do
seq=$(echo $line | cut -d' ' -f2)
cell=$(echo $line | cut -d' ' -f1)
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq/$cell
  srr_line=${line##"$cell $seq "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo "$srr $cell $seq"
    #star_rna
    #mark_duplicates
  done
else
  mark="$seq"
  if [[ "$seq" = GATA3* ]] || [[ "$seq" = PU.1 ]] || [[ "$seq" = Input ]]; then
    cd $dir0/Data/ChIP_seq/Transcription_Factor/$cell/$mark
    srr_line=${line##"$cell $seq "}
    IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo "$srr $cell $seq"
    #star_chip
    #mark_duplicates
  done
  else
    cd $dir0/Data/ChIP_seq/$cell/$mark
    srr_line=${line##"$cell $seq "}
    IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo "$srr $cell $seq"
    #star_chip
    #mark_duplicates
  done
  fi
fi
done < $infodir/srr_files_rothenberg.txt

if [ -f $paralleldir/mark_duplicates_rothenberg.txt ]; then
parallel --xapply --dryrun -j 40 -- < $paralleldir/mark_duplicates_rothenberg.txt
#parallel --xapply -j 40 -- < $paralleldir/mark_duplicates_rothenberg.txt
fi
if [ -f $paralleldir/remove_duplicates_rothenberg.txt ]; then
parallel --xapply --dryrun -j 40 -- < $paralleldir/remove_duplicates_rothenberg.txt
#parallel --xapply -j 40 -- < $paralleldir/remove_duplicates_rothenberg.txt
fi

##---------------------COMBINE BAM FILES ROTHENBERG---------------------------------------
#Combine fastq files, compress to .gz, and align
lines=$(wc -l $infodir/sample_files_rothenberg.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_rothenberg.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
seq=$(echo $line | cut -d':' -f2 | xargs | cut -d' ' -f2 | xargs)
cell=$(echo $line | cut -d':' -f2 | xargs | cut -d' ' -f1 | xargs)
rep=$(echo $line | cut -d':' -f2 | xargs | cut -d' ' -f3 | xargs)
echo "$seq"
if [[ "$seq" = RNA_seq* ]]; then
  cd $dir0/Data/RNA_seq/$cell
  filename=${cell}_${seq}_Rothenberg_mm10_${gsm}
  echo "$filename"
  cat SRR*_trimmed.fq > ${filename}.fastq
  pigz -f -p 40 ${filename}.fastq
  #heatmap
  bam_combine
else
  mark="$seq"
  if [[ "$seq" = GATA3* ]] || [[ "$seq" = PU.1 ]] || [[ "$seq" = Input ]]; then
    cd $dir0/Data/ChIP_seq/Transcription_Factor/$cell/$seq
    filename=${cell}_ChIP_seq_${seq}_Rothenberg_mm10_${gsm}
    echo "$filename"
    cat SRR*_trimmed.fq > ${filename}.fastq
    pigz -f -p 40 ${filename}.fastq
    #heatmap
    bam_combine
  else
    cd $dir0/Data/ChIP_seq/$cell/$mark
    filename=${cell}_ChIP_seq_${seq}_Rothenberg_mm10_${gsm}
    echo "$filename"
    cat SRR*_trimmed.fq > ${filename}.fastq
    pigz -f -p 40 ${filename}.fastq
    #heatmap
    bam_combine
  fi
fi
done

parallel --xapply --dryrun -j 30 -- < $paralleldir/merge_bam_rothenberg.txt
parallel --xapply -j 30 -- < $paralleldir/merge_bam_rothenberg.txt
#parallel --xapply --dryrun -j 30 -- < $paralleldir/heat_map_rothenberg.txt
#parallel --xapply -j 30 -- < $paralleldir/heat_map_rothenberg.txt
