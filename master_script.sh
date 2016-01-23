#!/bin/bash
start=`date +%s`
#rm -f log.log
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
mkdir -p $dir0/Analysis/ATAC_seq/R_output
mkdir -p $dir0/Analysis/ChIP_seq/R_output
mkdir -p $dir0/Analysis/MNase_seq/DP
mkdir -p $dir0/Analysis/Motif/Homer/Common_CD4_CD8_NK
mkdir -p $dir0/Analysis/ATAC_seq/Overlap
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_Discover/Cluster_Plots
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_FC1
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_FC75
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_Common/Intersect/Motif_Overlap
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_Common/NK_CD4_CD8
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_Common/NK_CD8
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_Common/Combine_All
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_Common/NK_CD4
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_Close/CD4_CD8
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_Close/NK_CD4_CD8
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_Intersect
mkdir -p $dir0/Analysis/ChIP_seq/Transcription_Factor/Motif/Tcf1/Homer_Motif_Tcf1
mkdir -p $dir0/Analysis/Homer/Tag_Directories/Amit
mkdir -p $dir0/Analysis/Homer/Tag_Directories/Dose
mkdir -p $dir0/Analysis/Homer/Tag_Directories/MNase
mkdir -p $dir0/Analysis/Homer/Tag_Directories/Rothenberg
mkdir -p $dir0/Data/MNase_seq/DP
mkdir -p $dir0/Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8_FC_NK
mkdir -p $datadir/MNase_seq/bam
mkdir -p $dir0/Data/MNase_seq/DP_Bfast


#--------------------Downloading Accessory Data-------------------------------------------
if [ ! -e $datadir/mm10.chrom.sizes ]; then
  cd $datadir && wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
  cd $dir0
fi

#---------------------Downloading Sample Info---------------------------------------------
#. ${scriptdir}/download_amit.sh
#. ${scriptdir}/download_tcf1.sh
#. ${scriptdir}/download_rothenberg.sh

#---------------------Downloading and Aligning Files to mm10------------------------------
#. ${scriptdir}/fastq-align-amit.sh

#---------------------Move Files and Generate MACS----------------------------------------
#. ${scriptdir}/move_and_macs_output.sh

##---------------------Normalize bigwig files from bam BIG WIG FILES FROM BAM-------------
#. ${scriptdir}/bam_coverage_bigwig.sh

#---------------------MACS Peak Calling---------------------------------------------------
#. ${scriptdir}/macs_peak_amit.sh

#---------------------Downloading Sample Info---------------------------------------------

#---------------------Downloading MNase DP------------------------------------------------
#. ${scriptdir}/download_mnase.sh

#---------------------Downloading and Aligning Files to mm10------------------------------
#. ${scriptdir}/fastq-align_tcf1.sh

#---------------------Downloading Sample Info---------------------------------------------

#---------------------Downloading and Aligning Files to mm10------------------------------
#. ${scriptdir}/fastq-align_rothenberg.sh

#---------------------Move Files and Generate MACS----------------------------------------
#. ${scriptdir}/move_and_macs_output_rothenberg.sh

#---------------------Normalization-------------------------------------------------------
#. ${scriptdir}/normalization.sh

#---------------------MACS Peaks Coverage-------------------------------------------------
#. ${scriptdir}/macs_coverage.sh

#---------------------Homer Motif TCF1----------------------------------------------------
#. ${scriptdir}/homer_tcf1_motif.sh

#---------------------Union Coverage Normalize------------------------------------------------------
#. ${scriptdir}/union_atac_all_normalize.sh

#---------------------Unique ATAC seq Peaks------------------------------------------------------
. ${scriptdir}/atac_seq_unique_all.sh

#---------------------MNase Macropahge------------------------------------------------------
. ${scriptdir}/MNase_Macrophage.sh

#---------------------Homer Motif Search Tcf1------------------------------------------------------
#. ${scriptdir}/homer_motif_search_tcf1.sh

#---------------------Union Coverage------------------------------------------------------
#. ${scriptdir}/atac_union_coverage.sh

#---------------------CD4 CD8 NK Common Motif Analysis------------------------------------
#. ${scriptdir}/motif_common_CD4_CD8_NK.sh

#---------------------Unique ATAC ChIP Coverage-------------------------------------------
#. ${scriptdir}/unique_atac_chip_coverage.sh
