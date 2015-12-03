#!/bin/bash
#This script will take Peak files from each ATAC seq sample (Amit, 2014) and perform bedtools union
#Tag coverage will then be performed for each cell at the coordinates from the union file
#Tag coverage will be normalized to number of reads per million

mkdir -p $dir0/Analysis/ATAC_seq/Union
cd $datadir/ATAC_seq/macs

#---------------------CREATE TEMP BED FILE------------------------------------------------
##This generates a temp bed file for unionbedg function later
#This temp file contains the ATAQ-seq peaks called by MACS (CHR, Start, End, and Value)
mkdir -p $dir0/Analysis/ATAC_seq/Union

cd $datadir/ATAC_seq/macs
for i in `ls *tag_peaks.bedgraph`;
do
  filename=`echo "$i" | cut -d'.' -f1`		 ##Generates filename variable
  cp $i ${dir0}/Analysis/ATAC_seq/Union/${filename}.bed		##Copies to temp file
done

#---------------------UNION FILE----------------------------------------------------------
#Creating combined peak file using files generated above using bedtools unionbedg
#--------------------OPTIONS--------------------------------------------------------------
cd $dir0/Analysis/ATAC_seq/Union
for i in `ls *tag_peaks.bed`; do
filename=`echo "$i" | cut -d'.' -f1`
cp $i ${filename}_temp.bed 		##Copies to temp file
done
filename=`ls *_temp.bed  | cut -d'_' -f1 | paste -sd "_" -` 	##Generates combined filename variable
echo "$filename"
bedtools unionbedg -i `ls *_temp.bed` > ${filename}_merged.bed ##Combines the peak files together

Rscript --verbose ${r_dir}/ATAC_Peaks_Union_Splitter.R $dir0 ##Splits combined peak file into separate peak files

#------------COVERAGE ATAC SEQ PEAKS------------------------------------------------------
#Gets tag coverage from bam files for each separate peak file generated above
rm -f $paralleldir/atac_coverage_parallel.txt
for i in `ls *union_split.bed`; do
  cell=${i%*_union_split.bed}
  bamfile=$(find $datadir/ATAC_seq/bam/${cell}*.bam)
  filename=$(echo "$bamfile" | cut -d'.' -f1 | rev | cut -d'/' -f1 | rev)
  echo `pwd`"/${filename}_Union_Split_Tag.bed" `pwd`"/${i}" ${bamfile} >> $paralleldir/atac_coverage_parallel.txt	
done

parallel --xapply --dryrun --colsep ' ' -a $paralleldir/atac_coverage_parallel.txt "bedtools coverage -sorted -g $datadir/sorted.mm10.genome -counts -a {2} -b {3} > {1}"
parallel --xapply --colsep ' ' -a $paralleldir/atac_coverage_parallel.txt "bedtools coverage -sorted -g $datadir/sorted.mm10.genome -counts -a {2} -b {3} > {1}"

#Removes old tag counts and keeps only the tag counts from bamtools coverage
for i in `ls *_Union_Split_Tag.bed`; 
do
    filename=${i%*.bed}
    awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5}' ${i} > ${filename}.bedgraph	##Removes the old tag counts, keeping only the new tag counts
done


#----------NORMALIZATION------------------------------------------------------------------
#Normalize each peak file
#Rscript --verbose ${r_dir}/Normalize_ATAC_Peaks_All.R $dir0	

cd $dir0/Analysis/ATAC_seq/Union
rm *_Union_Split_Tag.bedgraph
rm *_union_split.bed
rm *_tag_peaks.bed
rm *_tag_peaks_temp.bed
rm *_Union_Split_Tag.bed