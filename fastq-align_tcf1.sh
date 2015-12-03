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

#STAR alignment for ChIP-seq
function star_chip {
STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate 
mv Aligned.sortedByCoord.out.bam STAR_${srr}.bam #Keeps the bam file that STAR generates
mv Log.final.out ${srr}_Log.final.out #Keeps a copy of the aligner log
STAR --runMode inputAlignmentsFromBAM --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --inputBAMfile STAR_${srr}.bam --bamRemoveDuplicatesType UniqueIdentical
mv Processed.out.bam dedup_${srr}.bam #Keeps a bam file with PCR duplicates removed via STAR
}

function mark_duplicates {
if [[ $(head -1 ${srr}.fastq | cut -c 1-4) = "@SRR" ]]; then
  #If file does not have header information regarding location of the read on the lane, will not look for optical duplicates
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/STAR_${srr}.bam OUTPUT="`pwd`"/${srr}.bam READ_NAME_REGEX=null REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_tcf1.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/STAR_${srr}.bam OUTPUT="`pwd`"/dupsMarked_${srr}.bam READ_NAME_REGEX=null METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_tcf1.txt

else
  #If file has proper header information regarding location of the read on the lane, will look for optical duplicates
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/STAR_${srr}.bam OUTPUT="`pwd`"/${srr}.bam REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_tcf1.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/STAR_${srr}.bam OUTPUT="`pwd`"/dupsMarked_${srr}.bam METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_tcf1.txt
fi
}


#STAR remove genome from memory
function star_remove {
STAR --genomeLoad Remove --genomeDir /mnt/data0/John/genome_GRCm38p4_M6
}

#Make sure genome is unloaded on first-run
cd $dir0
star_remove
#rm -f $infodir/pcr_duplicates.txt
rm -f $paralleldir/mark_duplicates_tcf1.txt
rm -f $paralleldir/remove_duplicates_tcf1.txt
rm -f $paralleldir/merge_bam_tcf1.txt
rm -f $paralleldir/heat_map_tcf1.txt
rm -f $paralleldir/homer_key_file.txt
rm -f $paralleldir/macs_parallel_tf.txt

##---------------------DOWNLOAD SRR TCF-1 DATA--------------------------------------------
#Download the SRR files for TCF-1 Data
if [ -f $paralleldir/srr_download_tcf1.txt ]; then
parallel --xapply --dryrun -j 30 --colsep ' ' -a $paralleldir/srr_download_tcf1.txt  "fastq-dump {3} {1} -O {2}"
#parallel --xapply -j 30 --colsep ' ' -a $paralleldir/srr_download_tcf1.txt  "fastq-dump {3} {1} -O {2}"
fi

##---------------------ALIGN TCF-1 DATA---------------------------------------------------
#Align fastq to bam for TCF-1 data
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
    echo $srr "align and mark duplicates"
    #star_chip
    #mark_duplicates
  done
fi
done < $infodir/srr_files_tcf1.txt

if [ -f $paralleldir/mark_duplicates_tcf1.txt ]; then
parallel --xapply --dryrun -j 40 -- < $paralleldir/mark_duplicates_tcf1.txt
#parallel --xapply -j 40 -- < $paralleldir/mark_duplicates_tcf1.txt
fi
if [ -f $paralleldir/remove_duplicates_tcf1.txt ]; then
parallel --xapply --dryrun -j 40 -- < $paralleldir/remove_duplicates_tcf1.txt
#parallel --xapply -j 40 -- < $paralleldir/remove_duplicates_tcf1.txt
fi

##---------------------COMBINE BAM FILES TCF-1--------------------------------------------
#Combine fastq files, compress to .gz, combine bam files
lines=$(wc -l $infodir/sample_files_tcf1.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_tcf1.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
seq=$(echo $line | cut -d':' -f2 | cut -d';' -f3 | xargs)
cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs)
#echo $gsm
if [[ "$seq" = ChIP* ]]; then
  mark=$(echo $line | cut -d':' -f2 | xargs | cut -d';' -f1 | cut -d' ' -f1 | xargs)
  cell=$(echo $line | cut -d':' -f2 | xargs | cut -d';' -f1 | cut -d' ' -f2 | xargs)
  cd $dir0/Data/ChIP_seq/Transcription_Factor/$cell/$mark
  filename=${cell}_${seq}_${mark}_Dose_mm10_${gsm}
  #echo "$filename"
  #cat SRR*.fastq > ${filename}.fastq
  #pigz -f -p 40  ${filename}.fastq
  #cp "$(ls SRR*.bam)" ${filename}.bam
  #samtools index ${filename}.bam
  #cp ${filename}.bam $datadir/ChIP_seq/bam
  #cp ${filename}.bam.bai $datadir/ChIP_seq/bam
  #cp ${filename}.fastq.gz $datadir/ChIP_seq/fastq
  echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
fi
done

cd $dir0/Analysis/Homer/Tag_Directories/Dose
#batchMakeTagDirectory.pl $paralleldir/homer_key_file.txt -cpu 40 -genome mm10 -format sam

unset input
unset treat
unset treatFrag
unset name
##---------------------COMBINE BAM FILES TCF-1--------------------------------------------
#Combine fastq files, compress to .gz, combine bam files
lines=$(wc -l $infodir/sample_files_tcf1.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_tcf1.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
seq=$(echo $line | cut -d':' -f2 | cut -d';' -f3 | xargs)
cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs)
if [[ "$seq" = ChIP* ]]; then
  mark=$(echo $line | cut -d':' -f2 | xargs | cut -d';' -f1 | cut -d' ' -f1 | xargs)
  cell=$(echo $line | cut -d':' -f2 | xargs | cut -d';' -f1 | cut -d' ' -f2 | xargs)
  cd $dir0/Data/ChIP_seq/Transcription_Factor/$cell/$mark
  filename=${cell}_${seq}_${mark}_Dose_mm10_${gsm}
  fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/Dose/$filename/tagInfo.txt" | cut -d"=" -f2 | xargs)
  if [[ "$mark" = Input* ]]; then
  input=$(echo `pwd`/"$filename")
  echo "$filename"
  fi 
  if [ "$mark" != "Input" ]; then
  name=$filename
  treat=$(echo `pwd`/"$filename")
  treatFrag=$fragLength
  echo "$filename"
  echo "Generating bigWig for $filename"
  #bamCoverage --bam ${filename}.bam --outFileName ${filename}.bw --outFileFormat bigwig --bamIndex ${filename}.bam.bai --normalizeUsingRPKM --binSize 10 --centerReads --binSize 10 --numberOfProcessors 38 --fragmentLength $fragLength
  cp ${filename}.bw $datadir/ChIP_seq/bw
  fi
fi
done
echo "${treat}.bam" "${name}_macs" "${input}.bam" "$treatFrag" >> $paralleldir/macs_parallel_tf.txt

cd $dir0/Analysis/ChIP_seq/Transcription_Factor/macs
MACSpvalue=1e-7
#parallel --dryrun --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_tf.txt "macs14 -t {1} -n {2} -c {3} --bw {4} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"
#parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_tf.txt "macs14 -t {1} -n {2} -c {3} --bw {4} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"

#---------------------CREATE BEDGRAPH FROM PEAKS------------------------------------------
##This generates a bedgraph of peaks using number of tags from MACS output using the generated xls file 
for i in `ls *macs_peaks.xls`; do
  filename=`echo "$i" | cut -d'.' -f1`
  $scriptdir/Macs2Bedgraph ${filename}.xls > ${filename}.bedgraph	##Custom script to append columns of CHR, Start and End to number of tags from xls file
  #sed -i '/_random/d' ${filename}.bedgraph
  bamname=${i%*_macs_peaks.xls}
  find $datadir/ChIP_seq/bam/${bamname}.bam
  sed -e 's/chr1	/chr1A	/' ${filename}.bedgraph  | sort -k1,1 -k2,2n | sed -e 's/chr1A/chr1/' > ${filename}_sort.bedgraph 
  bedtools coverage -sorted -g $datadir/sorted.mm10.genome -counts -a ${filename}_sort.bedgraph  -b $(find $datadir/ChIP_seq/bam/${bamname}.bam) | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' - > ${filename}_tag_peaks.bedgraph
done

for i in `ls *.bedgraph`;
do
  filename=`echo "$i" | cut -d'.' -f1`		 ##Generates filename variable
  cp $i ${dir0}/Analysis/ATAC_seq/temp/${filename}.bed 		##Copies to temp file
  #--------------------Sorting file---------------------------------------------------------
  sed -e 's/chr1	/chr1A	/' ${dir0}/Analysis/ATAC_seq/temp/${filename}.bed  | sort -k1,1 -k2,2n | sed -e 's/chr1A/chr1/' > ${dir0}/Analysis/ATAC_seq/temp/${filename}_sort.bed 
done

#---------------------GENERATE MACS MODEL BUILDING FIGURE---------------------------------
#for i in `ls *model.r`; do
#  Rscript --verbose $i
#done

#---------------------Normalize WIG FILE--------------------------------------------------
#wd=`pwd`
#. ${scriptdir}/normalize_wig.sh

##----------------------Generate BigWig file ChIP seq-------------------------------------	
parallel --xapply --dryrun -j 38 --colsep ' ' -a $paralleldir/macs_parallel_tf.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig $datadir/mm10.chrom.sizes {2}.bw"	##Converts wig file to BigWig
#parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_tf.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig $datadir/mm10.chrom.sizes {2}.bw"
parallel --xapply --dryrun -j 38 --colsep ' ' -a $paralleldir/macs_parallel_tf.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_normalized.wig $datadir/mm10.chrom.sizes {2}_normalized.bw"	##Converts wig file to BigWig
#parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_tf.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_normalized.wig $datadir/mm10.chrom.sizes {2}_normalized.bw"

