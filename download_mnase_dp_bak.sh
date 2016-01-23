#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject

##---------------------SET BIOPROJECT NUMBER----------------------------------------------
#query=PRJNA150959

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

