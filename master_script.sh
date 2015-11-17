#!/bin/bash
start=`date +%s`
#rm -f *.log
#logfile=log.log
#exec > $logfile 2>&1

#---------------------Setting Directories-------------------------------------------------
maindir=/mnt/data1/John/Pioneer_Factors
scriptdir=/mnt/data1/John/Pioneer-Factors
filedir=/mnt/data1/John/Pioneer-Factors/sample_files
paralleldir=/mnt/data1/John/Pioneer-Factors/parallel_commands
r_dir=/mnt/data1/John/Pioneer-Factors/R_Scripts
infodir=/mnt/data1/John/Pioneer-Factors/info
datadir=/mnt/data1/VahediLab/PTF_Team/Data
dir0=$maindir
cd ${dir0}
date > ${dir0}/Time.txt
echo "Start time: $start" >> ${dir0}/Time.txt
#Put all bam files for the cells to be analyzed into a folder named Data

#---------------------Making Directories--------------------------------------------------
#Makes the directories
mkdir -p $dir0/Data/ATAC_seq
mkdir -p $dir0/Data/ChIP_seq
mkdir -p $dir0/Data/RNA_seq
mkdir -p $dir0/Data/ChIP_seq/Transcription_Factor/
#mkdir -p $datadir/ATAC_seq/fastq ; mkdir -p $datadir/ATAC_seq/bam ; mkdir -p $datadir/ATAC_seq/bw ; mkdir -p $datadir/ATAC_seq/macs ; mkdir -p $datadir/ChIP_seq/fastq ; mkdir -p $datadir/ChIP_seq/bam ; mkdir -p $datadir/ChIP_seq/bw ; mkdir -p $datadir/ChIP_seq/macs ; mkdir -p $datadir/RNA_seq/fastq ; mkdir -p $datadir/RNA_seq/bam ; mkdir -p $datadir/RNA_seq/bw
mkdir -p $dir0/Analysis/ATAC_seq/macs
mkdir -p $dir0/Analysis/ATAC_seq/srr_macs
mkdir -p $dir0/Analysis/ChIP_seq/macs
mkdir -p $dir0/Analysis/ChIP_seq/srr_macs
#mkdir -p $dir0/Analysis/ChIP_seq/macs2
mkdir -p $dir0/Analysis/ChIP_seq/Transcription_Factor/macs
mkdir -p ${dir0}/Analysis/ATAC_seq/R_output
mkdir -p ${dir0}/Analysis/ChIP_seq/R_output
mkdir -p ${dir0}/Analysis/Motif/Homer/Common_CD4_CD8_NK
mkdir -p $dir0/Analysis/Homer/Tag_Directories/Amit

#--------------------Downloading Accessory Data-------------------------------------------
if [ ! -e $datadir/mm10.chrom.sizes ]; then
  cd $datadir && wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
  cd $dir0
fi

#---------------------Downloading Sample Info---------------------------------------------
#. ${scriptdir}/download.sh

#---------------------Downloading and Aligning Files to mm10------------------------------
#. ${scriptdir}/fastq-align.sh

#---------------------Move Files and Generate MACS----------------------------------------
#. ${scriptdir}/move_and_macs_output.sh

##---------------------Normalize bigwig files from bam BIG WIG FILES FROM BAM-------------
#. ${scriptdir}/bam_coverage_bigwig.sh

#---------------------MACS Peak Calling---------------------------------------------------
. ${scriptdir}/macs_peak_amit.sh


#---------------------Normalization-------------------------------------------------------
#. ${scriptdir}/normalization.sh

#---------------------Union Coverage------------------------------------------------------
#. ${scriptdir}/atac_union_coverage.sh

#---------------------CD4 CD8 NK Common Motif Analysis------------------------------------
#. ${scriptdir}/motif_common_CD4_CD8_NK.sh

#---------------------Unique ATAC ChIP Coverage-------------------------------------------
#. ${scriptdir}/unique_atac_chip_coverage.sh
