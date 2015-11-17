#!/bin/bash
#This script will add 500bp on each side of the Tcf1 binding peak and get tag coverage for the histone modifications

#---------------------Setting Directories-------------------------------------------------
dir0=/mnt/data1/John/pioneer
scriptdir=/mnt/data1/John/Pioneer-Factors
filedir=/mnt/data1/John/Pioneer-Factors/files/
mkdir -p $dir0/analysis/chip_seq/tcf1/tag
mkdir -p $dir0/analysis/chip_seq/tcf1/macs/motif

#---------------------Setting Window around each peak-------------------------------------
##This generates a bedgraph of peaks with a 500bp window around the start and end
cd ${dir0}/analysis/chip_seq/tcf1/macs
for i in `ls  *macs_peaks.bedgraph`; do
  filename=`echo "$i" | cut -d'_' -f1,2`
  awk 'BEGIN{OFS="\t"} {start=$2-1000;end=$3+1000;print $1,start,end}' $i > ${dir0}/analysis/chip_seq/tcf1/tag/${filename}_window.bed
  awk 'BEGIN{OFS="\t"} {start=$2-1000;end=$3+500;print $1,start,end}' $i > ${dir0}/analysis/chip_seq/tcf1/tag/${filename}_window_500.bed
done
  
#--------------------Sorting file---------------------------------------------------------
sed -e 's/chr1	/chr1A	/' ${dir0}/analysis/chip_seq/tcf1/tag/${filename}_window.bed | sort -k1,1 -k2,2n | sed -e 's/chr1A/chr1/' > ${dir0}/analysis/chip_seq/tcf1/tag/${filename}_window_sort.bed	

#---------------------CREATE 6 COLUMN BED FILE--------------------------------------------
##This generates a 6 column bedgraph that includes unique peak names and adds strand info used for homer

awk 'BEGIN {OFS="\t"} {print $0 "\tPeak" NR "\t0" "\t+"}' ${filename}_window_sort.bed > ${filename}_window_sort_homer.bed
annotatePeaks.pl ${filename}_window_sort_homer.bed mm10 | sort -t $'\t' -k2,1 -k3,2n > ${filename}_window_sort_homer_annotate.bed
cat ${filename}_window_sort_homer_annotate.bed | sed -e 's/Peak//1' | sort -k1,1n > ${filename}_window_sort_homer_annotate_sort.bed

#------------COVERAGE HISTONE MODIFICATIONS TCF1------------------------------------------
cell=${dir0}/analysis/chip_seq/tag/tcf1_Thy_window_sort.bed
for i in `ls $dir0/data/chip_seq/bam/{B*.bam,CD4*.bam,CD8*.bam,NK*.bam,LTHSC*.bam,STHSC*.bam,MPP*.bam,CLP*.bam}`; do
  bamfile=`echo "$i" | rev | cut -d'/' -f1 | rev | cut -d'.' -f1`
  ##Gets histone mark tag counts from bam files  $cell is the cell of interest and is option -a $i is the bam file to be probed and is option -b
  ##Stores in a file named with destination file, then cell (with peaks) listed first, then bam file that was probed
  echo "${dir0}/analysis/chip_seq/tcf1/tag/${bamfile}_Tcf1_Thy_tag.bed ${dir0}/analysis/chip_seq/tcf1/tag/${filename}_window_sort.bed $i"		
done > $filedir/chip_coverage_parallel.txt

#Getting total number of tags for ChIP-seq peaks used for normalization later
cd ${dir0}/data/chip_seq/
for i in `ls -d {LTHSC/*,STHSC/*,MPP/*,CLP/*,B/*,CD4/*,CD8/*,NK/*}`; do
cd $i
cell=$(echo $i | cut -d'/' -f1)
mark=$(echo $i | cut -d'/' -f2)
echo -e $cell $mark $(awk 'BEGIN { FS = "|\t" } {if (NR == 9) {printf "%.5f\n", $2/1000000}}' Log.final.out)
cd ${dir0}/data/chip_seq/
done > $filedir/mapped_reads_histone.txt

parallel --xapply --dryrun --colsep ' ' -a $filedir/chip_coverage_parallel.txt "bedtools coverage -sorted -g ${dir0}/data/sorted.mm10.genome -counts -a {2} -b {3} > {1}"
parallel --xapply --colsep ' ' -a $filedir/chip_coverage_parallel.txt "bedtools coverage -sorted -g ${dir0}/data/sorted.mm10.genome -counts -a {2} -b {3} > {1}"

#	#Getting total tag counts for all ChIP seq histone modifications
#	rm -f ${dir0}/Analysis/ChIP_seq/read_mark.txt		##Removes any old read_mark.txt file
#	for i in `ls ${maindir}/Data/ChIP_seq/*.bam`;
#do
#		filename=`echo "$i" | cut -d'.' -f1`		##Gets the filename
#		mark=`echo "$filename" | cut -d'_' -f1`		##Gets the histone mark
#		cell=`echo "$filename" | cut -d'_' -f2,3`			##Gets the cell name
#		reads=`samtools idxstats $i | cut -d$'\t' -f3 | paste -sd+ | bc`		##Gets the total number of tags from the bamfile
#		#paste <(./progA) <(./progB)
#		paste <(echo -e "$mark\t$cell\t$reads") <(echo -e "$reads" | awk '{printf "%.5f\n", $1/1000000}')		##Stores it in a text file
