#!/bin/bash
#This script will make tag directories for homer used in Pioneer Factors

#---------------------Setting Directories-------------------------------------------------
maindir=/mnt/data1/John/Pioneer_Factors
dir0=/mnt/data0/John/Pioneer_Factors/homer
cd ${dir0}
	
#----------------------HOMER MAKE TAG DIRECTORY-------------------------------------------
##Creates tag directories file used for homer (ChIP-seq)
rm -f ${dir0}/ChIP_seq/ChIP_key_file.txt		##Removes any old ChIP_key_file
for i in `ls ${maindir}/Data/ChIP_seq/{*_B.bam,*_CD4.bam,*_CD8.bam,*_NK.bam,*_HSC.bam,*_MPP.bam,*_CLP.bam}`;		##Lists all ChIP seq bam files
do
	filename=`echo "$i" | rev | cut -d'/' -f1 | rev | cut -d'.' -f1`
	echo -e ${dir0}/ChIP_seq/${filename}'\t'${i}		##Creates keyfile for tag directories file for homer
done >> ${dir0}/ChIP_seq/ChIP_key_file.txt
##Creates tag directories (ChIP-seq)
batchMakeTagDirectory.pl ${dir0}/ChIP_seq/ChIP_key_file.txt -cpu 35

##Creates tag directories file used for homer (ATAC-seq)
rm -f ${dir0}/ATAC_seq/ATAC_key_file.txt		##Removes any old ChIP_key_file
mkdir -p ${maindir}/Analysis/ATAC_seq/homer
for i in `ls ${maindir}/Data/ATAC_seq/bam/*.bam`;		##Lists all ATAC seq bam files
do
	filename=`echo "$i" | rev | cut -d'/' -f1 | rev | cut -d'.' -f1`
	echo -e ${dir0}/ATAC_seq/${filename}'\t'${i}	##Creates tag directories file for homer
done >> ${dir0}/ATAC_seq/ATAC_key_file.txt
##Creates tag directories (ChIP-seq)
batchMakeTagDirectory.pl ${dir0}/ATAC_seq/ATAC_key_file.txt -cpu 35

#----------------------MEDIAN TAG COUNT---------------------------------------------------
##Will generate a file that contains the median tag count and average tag count for all ChIP seq bam fiels
cd ${dir0}/ChIP_seq
for i in `ls -d ${dir0}/ChIP_seq/H3*`;		##Lists all ChIP seq bam files
do
	filename=`echo "$i" | rev | cut -d'/' -f1 | rev`
	mark=`echo "$filename" | cut -d'_' -f1`	##Gets the histone mark
	cell=`echo "$filename" | cut -d'_' -f2,3`	##Gets the cell name
	median=$(head -n 1 ${dir0}/ChIP_seq/$filename/tagCountDistribution.txt | cut -d' ' -f7 | cut -d',' -f1)		##Gets the median tag count per bp
	average=$(sed -n '6p' ${dir0}/ChIP_seq/$filename/tagInfo.txt | cut -d$'\t' -f1 | cut -d'=' -f2)				##Gets the average tag count per bp
	echo -e "$cell\t$mark\t$median\t$average"		##Prints cell, histone mark, median, and average as one line
done > ChIP_median_tag_count.txt		##Prints all to a file
sort ChIP_median_tag_count.txt -o ChIP_median_tag_count.txt		##Sorts the file according to cell name

##Will generate a file that contains the median tag count and average tag count for all ATAC seq bam files
cd ${dir0}/ATAC_seq
for i in `ls -d ${dir0}/ATAC_seq/*/`;		##Lists all ATAC seq bam files
do
	filename=`echo "$i" | rev | cut -d'/' -f2| rev`
	median=$(head -n 1 ${dir0}/ATAC_seq/$filename/tagCountDistribution.txt | cut -d' ' -f7 | cut -d',' -f1)		##Gets the median tag count per bp
	average=$(sed -n '6p' ${dir0}/ATAC_seq/$filename/tagInfo.txt | cut -d$'\t' -f1 | cut -d'=' -f2)				##Gets the average tag count per bp
	echo -e "$filename\t$median\t$average"		##Prints cell, median, and average as one line
done > ATAC_median_tag_count.txt		##Prints all to a file
sort ATAC_median_tag_count.txt -o ATAC_median_tag_count.txt		##Sorts the file according to cell name

mkdir -p ${dir0}/R_output		##Creates R_output directories
Rscript --verbose ${maindir}/R_Scripts/MedianTag.R $dir0