#!/bin/bash
#This script will take ATAC_seq bam files, do MACS peak calling, create bedgraphs of peaks (based on # of tags), and create a file combining all peaks
#This script will take ChIP_seq bam files, do MACS peak calling, create bedgraphs of peaks (based on $ of tags), and create a file combining all peaks
#

#######---------------ATAC_SEQ COMBINED BAM-----------------------------------------------
#
#---------------------MACS PEAK CALLING---------------------------------------------------
##This uses MACS to determine ATAC_seq peaks on the replicate combined bam files
#-----------------------------------------------------------------------------------------
MACSpvalue=1e-7
cd ${dir0}/Analysis/ATAC_seq/macs
	
##----------------------MACS Peak Calling ATAC seq----------------------------------------
#parallel --xapply --dryrun -j 38 --colsep ' ' -a $paralleldir/macs_parallel_atac.txt "macs14 -t {1} -n {2} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"
#parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_atac.txt "macs14 -t {1} -n {2} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"

#---------------------GENERATE MACS MODEL BUILDING FIGURE---------------------------------
#for i in `ls *model.r`; do
#  Rscript --verbose $i
#done

#---------------------CHECK MACS ERRORS---------------------------------------------------
#wd=`pwd`
#. ${scriptdir}/macs_warning_test.sh

#---------------------Normalize WIG FILE--------------------------------------------------
#wd=`pwd`
#. ${scriptdir}/normalize_wig.sh



##----------------------Generate BigWig file ATAC seq-------------------------------------	
#parallel --xapply --dryrun -j 38 --colsep ' ' -a $paralleldir/macs_parallel_atac.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig $datadir/mm10.chrom.sizes {2}.bw"	##Converts wig file to BigWig
#parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_atac.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig $datadir/mm10.chrom.sizes {2}.bw" #&& rm -r *MACS_wiggle #&& rm *negative_peaks*
#parallel --xapply --dryrun -j 38 --colsep ' ' -a $paralleldir/macs_parallel_atac.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_normalized.wig $datadir/mm10.chrom.sizes {2}_normalized.bw"	##Converts wig file to BigWig
#parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_atac.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_normalized.wig $datadir/mm10.chrom.sizes {2}_normalized.bw"

#---------------------CREATE BEDGRAPH FROM PEAKS------------------------------------------
##This generates a bedgraph of peaks using number of tags from MACS output using the generated xls file 
#for i in `ls *.xls`; do
#  filename=`echo "$i" | cut -d'.' -f1`
#  $scriptdir/Macs2Bedgraph ${filename}.xls > ${filename}.bedgraph	##Custom script to append columns of CHR, Start and End to number of tags from xls file
#  #sed -i '/_random/d' ${filename}.bedgraph
#done
cp *.bedgraph $datadir/ATAC_seq/macs
#######---------------ATAC_SEQ BAM END----------------------------------------------------

#######---------------ChIP_SEQ COMBINED BAM AMIT------------------------------------------
#
#---------------------MACS PEAK CALLING---------------------------------------------------
##This uses MACS to determine ChIP_seq peaks on the replicate combined bam files
##----------------------ChIP seq Histone Mark Peak Calling--------------------------------
cd ${dir0}/Analysis/ChIP_seq/macs
#parallel --xapply --dryrun -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip.txt  "macs14 -t {1} -n {2} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"
#parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip.txt  "macs14 -t {1} -n {2} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"

#---------------------GENERATE MACS MODEL BUILDING FIGURE---------------------------------
#for i in `ls *model.r`; do
#  Rscript --verbose $i
#done
#---------------------CHECK MACS ERRORS---------------------------------------------------
#wd=`pwd`
#. ${scriptdir}/macs_warning_test.sh

#---------------------Normalize WIG FILE--------------------------------------------------
#wd=`pwd`
#. ${scriptdir}/normalize_wig.sh

##----------------------Generate BigWig file ChIP  seq-------------------------------------	
#parallel --xapply --dryrun -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig $datadir/mm10.chrom.sizes {2}_afterfiting_all.bw"	#Converts wig file to BigWig
#parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig $datadir/mm10.chrom.sizes {2}_afterfiting_all.bw"
#parallel --xapply --dryrun -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_normalized.wig $datadir/mm10.chrom.sizes {2}_normalized.bw"	#Converts wig file to BigWig
#parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_normalized.wig $datadir/mm10.chrom.sizes {2}_normalized.bw"
#
#---------------------CREATE BEDGRAPH FROM PEAKS------------------------------------------
##This generates a bedgraph of peaks using number of tags from MACS output using the generated xls file 
#for i in `ls *.xls`; do
#  filename=`echo "$i" | cut -d'.' -f1`
#  $scriptdir/Macs2Bedgraph ${filename}.xls > ${filename}.bedgraph	##Custom script to append columns of CHR, Start and End to number of tags from xls file
  #sed -i '/_random/d' ${filename}.bedgraph
#done
cp *.bedgraph $datadir/ChIP_seq/macs
#######---------------ChIP_SEQ AMIT END---------------------------------------------------

