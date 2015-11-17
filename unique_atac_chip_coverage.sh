#!/bin/bash
#----------------------------ChIP seq analysis--------------------------------------------
#This section of the script will take unique ATAC seq peaks and probe for histone modification data as well as doing motif analysis


#------------------Sorting BED Files to BAM File Order------------------------------------
cd ${dir0}/Analysis/ATAC_seq/R_output
for i in `ls  *_unique_all.bedgraph`;
do
  filename=`echo "$i" | cut -d'.' -f1`
  sed -e 's/chr1	/chr1A	/' ${filename}.bedgraph | sort -k1,1 -k2,2n | sed -e 's/chr1A/chr1/' > ${filename}_sorted.bedgraph ##Sorts file
done

#---------------------Setting Window around each peak-------------------------------------
##This generates a bedgraph of peaks with a 500bp window around the start and end
n=1
for i in `ls  *_unique_all_sorted.bedgraph`;
do
name[$((n++))]=`echo "$i" | cut -d'_' -f1`		##Stores names of cells into an array
filename=`echo "$i" | cut -d'_' -f1`
awk 'BEGIN{OFS="\t"} {start=$2-1000;end=$3+1000;print $1,start,end}' $i > ${dir0}/Analysis/ChIP_seq/${filename}_window.bed		##Opens the peaks by 500bp on each side
mkdir -p ${dir0}/Analysis/ChIP_seq/$filename		##Creates directories for later use
done
n=$((n-1))

#------------COVERAGE HISTONE MODIFICATIONS-----------------------------------------------
while read line; do	##Lists the ATAC seq bam files to do peaking calling  
  cell=$line
  for c in `ls ${dir0}/Data/mm10/ChIP_seq/bam/{B*.bam,CD4*.bam,CD8*.bam,NK*.bam,LTHSC*.bam,STHSC*.bam,MPP*.bam,CLP*.bam}`;
    do
      bamfile=`echo "$c" | rev | cut -d'/' -f1 | rev | cut -d'.' -f1`
      ##Gets histone mark tag counts from bam files
      ##$cell is the cell of interest and is option -a
      ##$c is the bam file to be probed and is option -b
      ##Stores in a file named with destination file, then cell (with peaks) listed first, then bam file that was probed
      echo "${dir0}/Analysis/ChIP_seq/$cell/${cell}_${bamfile}_bam_count.bedgraph ${dir0}/Analysis/ChIP_seq/${cell}_window.bed $c"		
    done > $filedir/${cell}_ChIP_coverage_parallel.txt
  parallel --xapply --dryrun --colsep ' ' -a $filedir/${cell}_ChIP_coverage_parallel.txt "bedtools coverage -sorted -g ${filedir}/mm10.genome -counts -a {2} -b {3} > {1}"
  parallel --xapply --colsep ' ' -a $filedir/${cell}_ChIP_coverage_parallel.txt "bedtools coverage -sorted -g ${filedir}/mm10.genome -counts -a {2} -b {3} > {1}"
done < ${filedir}/atac_cell.txt