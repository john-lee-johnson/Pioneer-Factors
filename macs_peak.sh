#!/bin/bash
#This script will take ATAC_seq bam files, do MACS peak calling, create bedgraphs of peaks (based on # of tags), and create a file combining all peaks
#This script will take ChIP_seq bam files, do MACS peak calling, create bedgraphs of peaks (based on $ of tags), and create a file combining all peaks
#This is the second script to be run

#---------------------Setting Directories-------------------------------------------------
dir0=/mnt/data1/John/Pioneer_Factors
scriptdir=/mnt/data1/John/Pioneer-Factors
paralleldir=/mnt/data1/John/Pioneer-Factors/parallel_commands


#######---------------ATAC_SEQ------------------------------------------------------------
#
#---------------------MACS PEAK CALLING---------------------------------------------------
##This uses MACS to determine ATAC_seq peaks in the provided bam files
#-----------------------------------------------------------------------------------------
MACSpvalue=1e-7
cd ${dir0}/Analysis/ATAC_seq/macs
	
##----------------------MACS Peak Calling ATAC seq----------------------------------------
parallel --xapply --dryrun --colsep ' ' -a $paralleldir/macs_parallel_atac.txt "macs14 -t {1} -n {2} -f BAM -g mm -p $MACSpvalue -w --single-profile"
#parallel --xapply --colsep ' ' -a $paralleldir/macs_parallel_atac.txt "macs14 -t {1} -n {2} -f BAM -g mm -p $MACSpvalue -w --single-profile"

##----------------------Generate BigWig file ATAC seq-------------------------------------	
parallel --xapply --dryrun --colsep ' ' -a $paralleldir/macs_parallel_atac.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz $datadir/mm10.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig
#parallel --xapply --colsep ' ' -a $paralleldir/macs_parallel_atac.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz $datadir/mm10.chrom.sizes {2}.bw" && rm -r *MACS_wiggle && rm *negative_peaks* && mv *.bw $datadir/ATAC_seq/bw
cp * $datadir/ATAC_seq/macs
#---------------------CREATE BEDGRAPH FROM PEAKS------------------------------------------
##This generates a bedgraph of peaks using number of tags from MACS output using the generated xls file 
for i in `ls *.xls`; do
  filename=`echo "$i" | cut -d'.' -f1`
  $scriptdir/Macs2Bedgraph ${filename}.xls > ${filename}.bedgraph	##Custom script to append columns of CHR, Start and End to number of tags from xls file
  #sed -i '/_random/d' ${filename}.bedgraph
done
#
#######---------------ATAC_SEQ END--------------------------------------------------------

#######---------------ChIP_SEQ AMIT-------------------------------------------------------
#
#---------------------MACS PEAK CALLING---------------------------------------------------
##This uses MACS to determine ChIP_seq peaks in the provided bam files
##----------------------ChIP seq Histone Mark Peak Calling--------------------------------
cd ${dir0}/Analysis/ChIP_seq/macs
for i in `ls $paralleldir/macs_parallel_chipa*`;do
parallel --xapply --dryrun --colsep ' ' -a $paralleldir/macs_parallel_chip.txt  "macs14 -t {1} -n {2} -f BAM -g mm -p $MACSpvalue -w --single-profile"
#parallel --xapply --dryrun --colsep ' ' -a $paralleldir/macs_parallel_chip.txt  "macs14 -t {1} -n {2} -f BAM -g mm -p $MACSpvalue -w --single-profile"
done
##----------------------Generate BigWig file ChIP  seq-------------------------------------	
parallel --xapply --dryrun --colsep ' ' -a $paralleldir/macs_parallel_chip.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz $datadir/mm10.chrom.sizes {2}_afterfiting_all.bw"	#Converts wig file to BigWig
#parallel --xapply --colsep ' ' -a $paralleldir/macs_parallel_chip.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz $datadir/mm10.chrom.sizes {2}_afterfiting_all.bw" && rm *negative_peaks* && mv *.bw $datadir/ChIP_seq/bw
#
#######---------------ChIP_SEQ AMIT END---------------------------------------------------

#######---------------ChIP_SEQ ROTHENBERG-------------------------------------------------
#
#---------------------MACS PEAK CALLING ROTHENBERG----------------------------------------
##This uses MACS to determine ChIP_seq peaks in the provided bam files
##----------------------ChIP seq Histone Mark Peak Calling--------------------------------
MACSpvalue=1e-4
parallel --xapply --dryrun --colsep ' ' -a $paralleldir/macs_parallel_chip_rothenberg.txt  "macs14 -t {1} -n {2} -f BAM -g mm -p $MACSpvalue -w --single-profile"
#parallel --xapply --colsep ' ' -a $paralleldir/macs_parallel_chip_rothenberg.txt  "macs14 -t {1} -n {2} -f BAM -g mm -p $MACSpvalue -w --single-profile"

##----------------------Generate BigWig file ChIP  seq-------------------------------------	
parallel --xapply --dryrun --colsep ' ' -a $paralleldir/macs_parallel_chip_rothenberg.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz $datadir/mm10.chrom.sizes {2}_afterfiting_all.bw"	#Converts wig file to BigWig
#parallel --xapply --colsep ' ' -a $paralleldir/macs_parallel_chip_rothenberg.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz $datadir/mm10.chrom.sizes {2}_afterfiting_all.bw" && rm -r *MACS_wiggle && rm *negative_peaks*
#mv *.bw $datadir/ChIP_seq/bw
#
#######---------------ChIP_SEQ ROTHENBERG END---------------------------------------------

#######----------------------TRANSCRIPTION FACTORS CHIP SEQ-------------------------------
##----------------------------------------------------------------------------------------

cd $dir0/Analysis/ChIP_seq/Transcription_Factor/macs
##----------------------MACS Peak Calling TF----------------------------------------------
parallel --xapply --dryrun --colsep ' ' -a $paralleldir/macs_parallel_tf.txt "macs14 -t {1} -c {3} -n {2} -f BAM -g mm -p $MACSpvalue -w --single-profile"
#parallel --xapply --colsep ' ' -a $paralleldir/macs_parallel_tf.txt "macs14 -t {1} -c {3} -n {2} -f BAM -g mm -p $MACSpvalue -w --single-profile"

##----------------------Generate BigWig file TF-------------------------------------------	
parallel --xapply --dryrun --colsep ' ' -a $paralleldir/macs_parallel_tf.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz $datadir/mm10.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig
#parallel --xapply --colsep ' ' -a $paralleldir/macs_parallel_tf.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz $datadir/mm10.chrom.sizes {2}_afterfiting_all.bw"	&& rm -r *MACS_wiggle && rm *negative_peaks*
#mv *.bw $datadir/ChIP_seq/bw
#-----------------------------------------------------------------------------------------

#---------------------CREATE BEDGRAPH FROM PEAKS------------------------------------------
##This generates a bedgraph of peaks using number of tags from MACS output using the generated xls file 
cd $dir0/Analysis/ChIP_seq/tf/macs 
for i in `ls *.xls`; do
  filename=`echo "$i" | rev | cut -d'.' -f2,3 | rev`
  $scriptdir/Macs2Bedgraph ${filename}.xls > ${filename}.bedgraph	##Custom script to append columns of CHR, Start and End to number of tags from xls file
  #sed -i '/_random/d' ${filename}.bedgraph
done
#
#######----------------------TRANSCRIPTION FACTORS CHIP SEQ END---------------------------

#cd ${dir0}/Analysis/ChIP_seq/macs
#cp * $datadir/ChIP_seq/macs
#cd ${dir0}/Analysis/ChIP_seq/Transcription_Factor/macs
