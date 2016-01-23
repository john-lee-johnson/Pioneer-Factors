#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory

echo "#--------------CARRYING OUT THE ACTUAL DOWNLOAD----------"
echo "#--------------DOWNLOADING-------------------------" >> $logdir/download_log_${des}.txt
parallel --xapply --dryrun -j 30 -- < $paralleldir/srr_download_${des}.txt >> $logdir/download_log_${des}.txt
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
parallel --xapply --dryrun -j 30 -- < $paralleldir/srr_download_${des}.txt
parallel --xapply -j 30 -- < $paralleldir/srr_download_${des}.txt
echo "-------------------------------------------------------"; echo ""
echo "******Recommended to stop here and check a few of the files for quality and adapter contamination using FASTQC******"; echo ""