#!/bin/bash
#This script will make tag directories for homer used in Pioneer Factors

#---------------------Setting Directories-------------------------------------------------
maindir=/mnt/data1/John/Pioneer_Factors
dir0=/mnt/data0/John/Pioneer_Factors/homer
cd ${dir0}
	
#----------------------HOMER MAKE TAG DIRECTORY-------------------------------------------
##Creates tag directories file used for homer (ChIP-seq)
rm -f ${dir0}/ChIP_seq/ChIP_key_file.txt		##Removes any old ChIP_key_file
for i in `ls ${maindir}/Data/ChIP_seq/H3*.bam`;		##Lists all ChIP seq bam files
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