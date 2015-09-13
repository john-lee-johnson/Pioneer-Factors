#!/bin/bash
#This script will take ATAC_seq bam files, do MACS peak calling, create bedgraphs of peaks (based on # of tags), and create a file combining all peaks

#---------------------Setting Directories-------------------------------------------------
dir0=/mnt/data1/John/pioneer
filedir=/mnt/data1/John/Pioneer-Factors/files/
mkdir -p $dir0/analysis/atac_seq/macs
mkdir -p $dir0/analysis/chip_seq/tcf1
#cd ${dir0}/data && wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

#---------------------MACS PEAK CALLING---------------------------------------------------
##This uses MACS to determine ATAC_seq peaks in the provided bam files
#--------------------OPTIONS--------------------------------------------------------------
option=1	#Set to 1 if you want to run MACS peak calling
#-----------------------------------------------------------------------------------------
MACSpvalue=1e-7
cd ${dir0}/analysis/atac_seq/macs
	
##----------------------MACS Peak Calling ATAC seq----------------------------------------
parallel --xapply --dryrun --colsep ' ' -a $filedir/macs_parallel_atac.txt "macs14 -t {1} -n {2} -f BAM -g mm -s {4} -p {3} -w --single-profile"
#parallel --xapply --colsep ' ' -a $filedir/macs_parallel_atac.txt "macs14 -t {1} -n {2} -f BAM -g mm -s {4} -p {3} -w --single-profile"

##----------------------Generate BigWig file ATAC seq-------------------------------------	
parallel --xapply --dryrun --colsep ' ' -a $filedir/macs_parallel_atac.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz ${dir0}/data/mm10.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig
parallel --xapply --colsep ' ' -a $filedir/macs_parallel_atac.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz ${dir0}/data/mm10.chrom.sizes {2}.bw"	##Converts wig file to BigWig

cd $dir0/analysis/chip_seq/tcf1
##----------------------MACS Peak Calling Tcf1--------------------------------------------
parallel --xapply --dryrun --colsep ' ' -a $filedir/macs_parallel_tcf1.txt "macs14 -t {1} -c {5} -n {2} -f BAM -g mm -s {4} -p {3} -w --single-profile"
#parallel --xapply --colsep ' ' -a $filedir/macs_parallel.txt "macs14 -t {1} -n {2} -p {3} -w --single-profile"

##----------------------Generate BigWig file Tcf1-----------------------------------------	
parallel --xapply --dryrun --colsep ' ' -a $filedir/macs_parallel_tcf1.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz ${dir0}/data/mm10.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig
#parallel --xapply --dryrun --colsep ' ' -a $filedir/macs_parallel_tcf1.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz ${dir0}/data/mm10.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig
	
##----------------------Generate BigWig file-----------------------------------------	
parallel --xapply --dryrun --colsep ' ' -a $filedir/macs_parallel_atac.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz ${dir0}/data/mm10.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig
#parallel --xapply --colsep ' ' -a ${dir0}/options/macs_parallel.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz ${dir0}/data/mm10.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig	
#parallel --xapply --colsep ' ' -a ${dir0}/options/macs_parallel.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz ${dir0}/Data/mm9.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig	
#	rm -r *_macsoutput_MACS_wiggle		##Removes MACS output wig directory	
#	#Getting total number of tags for ATAC-seq peaks used for normalization later
#	cd ${dir0}/Analysis/ATAC_seq
#	rm -f ${dir0}/Analysis/ATAC_seq/read_count.txt
#	for i in `ls *.xls`;			##Lists MACS output xlx file
#	do
#		filename=`echo "$i" | cut -d'.' -f1`
#		reads=`awk 'BEGIN { FS = ":" } {if (NR == 15) {printf "%.5f\n", $2/1000000}}' ${filename}.xls`			##Gets total number of tags and divides by 10^6 and stores it in an array for later
#		cell=`echo $filename | cut -d'_' -f1`
#		echo -e "$cell\t$reads"
#	done >> ${dir0}/Analysis/ATAC_seq/read_count.txt