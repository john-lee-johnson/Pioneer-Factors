#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject

##---------------------SET BIOPROJECT NUMBER----------------------------------------------
query=PRJNA243033

##---------------------SET SAMPLE DESCRIPTION---------------------------------------------
#Set description string that will identify the info and parallel command files being output in this script
#Examples: mnase_bmdm, amit, mnase_dp, etc
#des=mnase_dp

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
#------------------------BFAST 0.7.0a-----------------------------------------------------

cd $dir0/Data/MNase_seq/DP_Bfast
#cp $dir0/Data/MNase_seq/DP/SRR1209650_F3.csfasta $dir0/Data/MNase_seq/DP_Bfast
#cp $dir0/Data/MNase_seq/DP/SRR1209650_F3_QV.qual $dir0/Data/MNase_seq/DP_Bfast
echo "#---BFast Convert Reads---"
#solid2fastq -n 10000000 -o reads *.csfasta *.qual
echo "#---BFast Convert the reference (nucleotide space)---"
#cp /mnt/data0/John/GRCm38p4_mm10_pa_only.fa $dir0/Data/MNase_seq/DP_Bfast/mm10.fa
#bfast fasta2brg -f mm10.fa
echo "#---BFast Convert the reference (color space)---"
#bfast fasta2brg -f mm10.fa -A 1
echo "#---BFast Create the indexes---"
bfast index -f mm10.fa -n 32 -m 1111111111111111111111 -w 14 -i 1 -A 1
#bfast index -f mm10.fa -n 38 -m 111110100111110011111111111 -w 14 -i 2 -A 1
#bfast index -f mm10.fa -n 38 -m 10111111011001100011111000111111 -w 14 -i 3 -A 1
#bfast index -f mm10.fa -n 38 -m 1111111100101111000001100011111011 -w 14 -i 4 -A 1
#bfast index -f mm10.fa -n 38 -m 111111110001111110011111111 -w 14 -i 5 -A 1
#bfast index -f mm10.fa -n 38 -m 11111011010011000011000110011111111 -w 14 -i 6 -A 1
#bfast index -f mm10.fa -n 38 -m 1111111111110011101111111 -w 14 -i 7 -A 1
#bfast index -f mm10.fa -n 38 -m 111011000011111111001111011111 -w 14 -i 8 -A 1
#bfast index -f mm10.fa -n 38 -m 1110110001011010011100101111101111 -w 14 -i 9 -A 1
#bfast index -f mm10.fa -n 38 -m 111111001000110001011100110001100011111 -w 14 -i 10 -A 1
echo "#---BFast Search the indexes---"
#bfast match -f mm10.fa -A 1 -r reads.1.fastq > bfast.matches.file.mm10.1.bmf


#-----------------------------------------------------------------------------------------
cd $dir0/Data/MNase_seq/DP
filename="${cell}_${seq}_mm10_${gsm}"
echo ${filename}.bam
#mv MNase_DP_mapped_dupsRemoved.bam ${filename}.bam
#samtools index ${filename}.bam
#cp ${filename}.bam $datadir/MNase_seq/bam
#cp ${filename}.bam.bai $datadir/MNase_seq/bam
#cp SRR1209650_F3.csfasta $datadir/MNase_seq/reads
#cp SRR1209650_F3_QV.qual $datadir/MNase_seq/reads
echo -e "${filename}"'\t'`pwd`"/${filename}_start_extend_sorted.bam" >> $paralleldir/homer_key_file.txt
cd $dir0/Analysis/Homer/Tag_Directories/MNase
#batchMakeTagDirectory.pl $paralleldir/homer_key_file.txt -cpu 40 -genome mm10 -format sam

cd $dir0/Data/MNase_seq/DP
fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/MNase/$filename/tagInfo.txt" | cut -d"=" -f2)
echo "Generating bigWig for $filename $fragLength"
#bamCoverage --bam ${filename}.bam --outFileName ${filename}.bw --outFileFormat bigwig --bamIndex ${filename}.bam.bai --normalizeUsingRPKM --binSize 10 --centerReads --binSize 10 --numberOfProcessors 38 --fragmentLength $fragLength
#cp ${filename}.bw $datadir/MNase_seq/bw

cd $dir0/Analysis/MNase_seq/DP
cp $datadir/ChIP_seq/macs/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph .
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph > Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed ##Adds a unique peak name and strand info

#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DP_MNase_Seq_mm10_GSM1359852 > histfile.txt
#-----------------------------------------------------------------------------------------
#
#                        MOUSE ESC MNase
#
#-----------------------------------------------------------------------------------------

##---------------------SET AMIT BIOPROJECT NUMBER-----------------------------------------
query=PRJNA175501


##---------------------SETTING FUNCTIONS--------------------------------------------------
#Creates parallel commands for fastq-dump
function fastq-download {
  fastq-dump -F -X 1 -Z ${srr} > $infodir/test.txt #Will download one line of the file to check for header information
    if [[ $(head $infodir/test.txt | sed -n 1p | cut -d":" -f1) = "@1" ]]; then
      echo "$srr" `pwd` >> $paralleldir/srr_download_mnase_esc.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "$srr" `pwd` "-F" >> $paralleldir/srr_download_mnase_esc.txt #If file does have proper header information, will keep the original header
    fi
#echo `pwd` `pwd`"/${srr}.fastq" >> $paralleldir/srr_download_mnase_esc.txt #Will remove sequencing adapters
}

rm -f $infodir/srr_download_mnase_esc.txt
rm -f $infodir/srr_files_mnase_esc.txt
rm -f $paralleldir/srr_download_mnase_esc.txt
rm -f $paralleldir/trim_galore_mnase_esc.txt
rm -f $paralleldir/homer_key_file.txt

##---------------------DOWNLOAD MNase DP BIOPROJECT SAMPLE DATA---------------------------
esearch -db sra -query $query | efetch --format docsum | xtract -pattern DocumentSummary -element Title > $infodir/sample_files_mnase_esc.txt
sed -i -e 's/MNase-Seq/MNase_Seq/g' $infodir/sample_files_mnase_esc.txt
sed -i -e 's/ merged mononucleosomes//g' $infodir/sample_files_mnase_esc.txt
sed '' $infodir/sample_files_mnase_esc.txt

##---------------------DOWNLOAD SRR TCF1 DATA---------------------------------------------
#Read in samples and generate a file of SRR numbers for TCF1 Data
lines=$(wc -l $infodir/sample_files_mnase_esc.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_mnase_esc.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
seq=$(echo $line | cut -d':' -f2 | cut -d';' -f3 | xargs)
cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs)
cd $dir0/Data/MNase_seq ; mkdir -p $cell ; cd $cell
echo "$gsm" "$seq" "$cell" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_mnase_esc.txt
done

##---------------------SRR TCF1 DATA PARALLEL OUTPUT--------------------------------------
#Read in a file and set up parallel commands for later SRR download
while read line; do
  gsm=$(echo $line | cut -d' ' -f1)
  seq=$(echo $line | cut -d' ' -f2)
  cell=$(echo $line | cut -d' ' -f3)
  cd $dir0/Data/MNase_seq/$cell
  srr_line=${line##"$gsm $seq $cell "}
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo "$cell $srr download"
    fastq-download
  done
done < $infodir/srr_files_mnase_esc.txt

##---------------------SETTING FUNCTIONS--------------------------------------------------

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
rm -f $paralleldir/mark_duplicates_tcf1.txt
rm -f $paralleldir/remove_duplicates_tcf1.txt
rm -f $paralleldir/merge_bam_tcf1.txt
rm -f $paralleldir/heat_map_tcf1.txt
rm -f $paralleldir/homer_key_file.txt
rm -f $paralleldir/macs_parallel_tf.txt

##---------------------DOWNLOAD SRR TCF-1 DATA--------------------------------------------
#Download the SRR files for TCF-1 Data
if [ -f $paralleldir/srr_download_mnase_esc.txt ]; then
parallel --xapply --dryrun -j 30 --colsep ' ' -a $paralleldir/srr_download_mnase_esc.txt  "fastq-dump {3} {1} -O {2}" > $paralleldir/srr_download_mnase_esc_commands.txt
cd /mnt/data0/ncbi/public/sra
#rm *.sra.cache
while read line; do
  echo "$line"
#  `$line`
  #rm *.sra
done < $paralleldir/srr_download_mnase_esc_commands.txt
#parallel --xapply -j 30 --colsep ' ' -a $paralleldir/srr_download_mnase_esc.txt  "fastq-dump {3} {1} -O {2}"
fi

##---------------------SET AMIT BIOPROJECT NUMBER-----------------------------------------
query=SAMN02712184

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
      abi-dump -F "$srr"
}

rm -f $infodir/srr_files_mnase_dp.txt
rm -f $paralleldir/srr_download_mnase_dp.txt
rm -f $paralleldir/trim_galore_mnase_dp.txt
rm -f $paralleldir/homer_key_file.txt

##---------------------DOWNLOAD MNase DP BIOPROJECT SAMPLE DATA---------------------------
esearch -db sra -query $query | efetch --format docsum | xtract -pattern DocumentSummary -element Title > $infodir/sample_files_mnase_dp.txt
sed -i -e 's/MNase-Seq/MNase_Seq/g' $infodir/sample_files_mnase_dp.txt

##---------------------DOWNLOAD MNase DATA___---------------------------------------------
#Read in samples and generate a file of SRR numbers for Amit Data
line=$(sed -n "1p" < $infodir/sample_files_mnase_dp.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
seq=$(echo $line | cut -d':' -f2 | cut -d';' -f3 | xargs)
cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs | cut -d' ' -f1 | xargs)
echo "$gsm $seq $cell"
echo "$gsm $seq $cell" `esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_mnase_dp.txt


line=$(sed -n "1p" < $infodir/srr_files_mnase_dp.txt)
gsm=$(echo $line | cut -d' ' -f1 | xargs)
seq=$(echo $line | cut -d' ' -f2 | xargs)
cell=$(echo $line | cut -d' ' -f3 | xargs)
srr=$(echo $line | cut -d' ' -f4 | xargs)
echo "$gsm $seq $cell $srr"
cd $dir0/Data/MNase_seq/DP
#fastq-download
#bowtie -p 38 -S -C mm10_colorspace -f SRR1209650_F3.csfasta -Q SRR1209650_F3_QV.qual MNase_DP.sam
#samtools view -b -F 4 MNase_DP.sam > MNase_DP_mapped.bam
#java -Xmx2g -jar $PICARD SortSam INPUT=MNase_DP_mapped.bam OUTPUT=MNase_DUP_sorted.bam SORT_ORDER=coordinate
#java -Xmx2g -jar $PICARD MarkDuplicates INPUT=MNase_DUP_sorted.bam OUTPUT=MNase_DP_mapped_dupsRemoved.bam READ_NAME_REGEX=null REMOVE_DUPLICATES=true METRICS_FILE=MNase_DP_mapped_metrics.txt
#java -Xmx2g -jar $PICARD MarkDuplicates INPUT=MNase_DUP_sorted.bam OUTPUT=MNase_DP_mapped_dupsMarked.bam READ_NAME_REGEX=null METRICS_FILE=MNase_DP_mapped_metrics.txt
#bedtools bamtobed -i DP_MNase_Seq_mm10_GSM1359852.bam > DP_MNase_Seq_mm10_GSM1359852.bed

#---------------------generating READ STARTS file-------------------------------
#TSS are found by taking the first nucleotide of the gene (start for + strand, and last for - strand)
#awk 'OFS="\t"{if ($6=="+") {start=$2; end=start+1} else if ($6=="-") {start=$3-1;end=$3} print $1,start,end,$4,$5,$6}' DP_MNase_Seq_mm10_GSM1359852.bed > DP_MNase_Seq_mm10_GSM1359852_start.bed
#---------------------generating READ STARTS file-------------------------------
#TSS are found by taking the first nucleotide of the gene (start for + strand, and last for - strand)
#awk 'OFS="\t"{if ($6=="+") {start=$2+50; end=$2+100} else if ($6=="-") {start=$3-100;end=$3-50} print $1,start,end,$4,$5,$6}' DP_MNase_Seq_mm10_GSM1359852_start.bed > DP_MNase_Seq_mm10_GSM1359852_start_extend.bed
#sort -k1,1 -k2,2n DP_MNase_Seq_mm10_GSM1359852_start_extend.bed > DP_MNase_Seq_mm10_GSM1359852_start_extend_sorted.bed
#bedToBam -i DP_MNase_Seq_mm10_GSM1359852_start_extend_sorted.bed -g $datadir/mm10.chrom.sizes > DP_MNase_Seq_mm10_GSM1359852_start_extend_sorted.bam



#-----------------------------------------------------------------------------------------
cd $dir0/Data/MNase_seq/DP
filename="${cell}_${seq}_mm10_${gsm}"
echo ${filename}.bam
#mv MNase_DP_mapped_dupsRemoved.bam ${filename}.bam
#samtools index ${filename}.bam
#cp ${filename}.bam $datadir/MNase_seq/bam
#cp ${filename}.bam.bai $datadir/MNase_seq/bam
#cp SRR1209650_F3.csfasta $datadir/MNase_seq/reads
#cp SRR1209650_F3_QV.qual $datadir/MNase_seq/reads
echo -e "${filename}"'\t'`pwd`"/${filename}_start_extend_sorted.bam" >> $paralleldir/homer_key_file.txt
cd $dir0/Analysis/Homer/Tag_Directories/MNase
#batchMakeTagDirectory.pl $paralleldir/homer_key_file.txt -cpu 40 -genome mm10 -format sam

cd $dir0/Data/MNase_seq/DP
fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/MNase/$filename/tagInfo.txt" | cut -d"=" -f2)
echo "Generating bigWig for $filename $fragLength"
#bamCoverage --bam ${filename}.bam --outFileName ${filename}.bw --outFileFormat bigwig --bamIndex ${filename}.bam.bai --normalizeUsingRPKM --binSize 10 --centerReads --binSize 10 --numberOfProcessors 38 --fragmentLength $fragLength
#cp ${filename}.bw $datadir/MNase_seq/bw

cd $dir0/Analysis/MNase_seq/DP
cp $datadir/ChIP_seq/macs/Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph .
awk 'BEGIN {OFS="\t"} {print $1,$2,$3 "\tpeak" NR "\t"$4 "\t+"}' Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_macs_peaks_tag_peaks.bedgraph > Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed ##Adds a unique peak name and strand info

#annotatePeaks.pl Thy_ChIP_seq_Tcf1_Dose_mm10_GSM1133644_5_column.bed mm10 -size 4000 -hist 10 -d $dir0/Analysis/Homer/Tag_Directories/MNase/DP_MNase_Seq_mm10_GSM1359852 > histfile.txt


