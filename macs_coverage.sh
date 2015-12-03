#!/bin/bash
mkdir -p ${dir0}/Analysis/ATAC_seq/temp
mkdir -p ${dir0}/Analysis/ChIP_seq/temp
#---------------------CREATE TEMP BED FILE------------------------------------------------
##This generates a temp bed file for tag coverage
cd $dir0/Analysis/ATAC_seq/macs
for i in `ls *.bedgraph`;
do
  filename=`echo "$i" | cut -d'.' -f1`		 ##Generates filename variable
  cp $i ${dir0}/Analysis/ATAC_seq/temp/${filename}.bed 		##Copies to temp file
  #--------------------Sorting file---------------------------------------------------------
  sed -e 's/chr1	/chr1A	/' ${dir0}/Analysis/ATAC_seq/temp/${filename}.bed  | sort -k1,1 -k2,2n | sed -e 's/chr1A/chr1/' > ${dir0}/Analysis/ATAC_seq/temp/${filename}_sort.bed 
done

#------------COVERAGE ATAC SEQ PEAKS------------------------------------------------------
cd $dir0/Analysis/ATAC_seq/temp
rm -f $paralleldir/atac_coverage_parallel.txt
for i in `ls *_sort.bed`; do
  filename=${i%*_macs_peaks_sort.bed}
  echo `pwd`"/${filename}_Tag.bed" `pwd`"/${i}" $(find $datadir/ATAC_seq/bam/${filename}.bam) >> $paralleldir/atac_coverage_parallel.txt
  echo "$filename"	
done

parallel --xapply --dryrun --colsep ' ' -a $paralleldir/atac_coverage_parallel.txt "bedtools coverage -sorted -g $datadir/sorted.mm10.genome -counts -a {2} -b {3} > {1}"
parallel --xapply --colsep ' ' -a $paralleldir/atac_coverage_parallel.txt "bedtools coverage -sorted -g $datadir/sorted.mm10.genome -counts -a {2} -b {3} > {1}"

cd $dir0/Analysis/ATAC_seq/temp
for i in `ls *Tag.bed`;
do
  filename=${i%*_Tag.bed}		 ##Generates filename variable
  awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' $i > $datadir/ATAC_seq/macs/${filename}_tag_peaks.bedgraph
done

#---------------------CREATE TEMP BED FILE------------------------------------------------
##This generates a temp bed file for tag coverage
cd $dir0/Analysis/ChIP_seq/macs
for i in `ls *.bedgraph`;
do
  filename=`echo "$i" | cut -d'.' -f1`		 ##Generates filename variable
  cp $i ${dir0}/Analysis/ChIP_seq/temp/${filename}.bed 		##Copies to temp file
  #--------------------Sorting file---------------------------------------------------------
  sed -e 's/chr1	/chr1A	/' ${dir0}/Analysis/ChIP_seq/temp/${filename}.bed  | sort -k1,1 -k2,2n | sed -e 's/chr1A/chr1/' > ${dir0}/Analysis/ChIP_seq/temp/${filename}_sort.bed 
done

#------------COVERAGE ChIP SEQ PEAKS------------------------------------------------------
cd $dir0/Analysis/ChIP_seq/temp
rm -f $paralleldir/chip_coverage_parallel.txt
for i in `ls *_sort.bed`; do
  filename=${i%*_macs_peaks_sort.bed}
  echo `pwd`"/${filename}_Tag.bed" `pwd`"/${i}" $(find $datadir/ChIP_seq/bam/${filename}.bam) >> $paralleldir/chip_coverage_parallel.txt
  echo "$filename"	
done

parallel --xapply --dryrun --colsep ' ' -a $paralleldir/chip_coverage_parallel.txt "bedtools coverage -sorted -g $datadir/sorted.mm10.genome -counts -a {2} -b {3} > {1}"
parallel --xapply --colsep ' ' -a $paralleldir/chip_coverage_parallel.txt "bedtools coverage -sorted -g $datadir/sorted.mm10.genome -counts -a {2} -b {3} > {1}"

cd $dir0/Analysis/ChIP_seq/temp
for i in `ls *Tag.bed`;
do
  filename=${i%*_Tag.bed}		 ##Generates filename variable
  awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' $i > $datadir/ChIP_seq/macs/${filename}_tag_peaks.bedgraph
done