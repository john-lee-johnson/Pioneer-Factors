#!/bin/bash
#######---------------ATAC_SEQ SRR--------------------------------------------------------
#
#---------------------MACS PEAK CALLING---------------------------------------------------
##This uses MACS to determine ATAC_seq peaks on the individual SRR bam files
#-----------------------------------------------------------------------------------------
MACSpvalue=1e-7
cd ${dir0}/Analysis/ATAC_seq/srr_macs	
##----------------------MACS Peak Calling SRR ATAC seq------------------------------------
for i in $(ls -d *); do
  cd ${dir0}/Analysis/ATAC_seq/srr_macs/$i
  #parallel --xapply --dryrun -j 38 --colsep ' ' -a $paralleldir/macs_parallel_atac_srr_${i}.txt "macs14 -t {1} -n {2} -f BAM -g mm -p $MACSpvalue -w --single-profile"
  #parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_atac_srr_${i}.txt "macs14 -t {1} -n {2} -f BAM -g mm -p $MACSpvalue -w --single-profile"
  #parallel --xapply --dryrun -j 38 --colsep ' ' -a $paralleldir/macs_parallel_atac_srr_${i}.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz $datadir/mm10.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig
  #parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_atac_srr_${i}.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz $datadir/mm10.chrom.sizes {2}.bw" #&& rm -r *MACS_wiggle && rm *negative_peaks*
  #---------------------GENERATE MACS MODEL BUILDING FIGURE---------------------------------
  #for a in `ls *model.r`; do
  #  Rscript --verbose $a
  #done
done
  
#######---------------ATAC_SEQ SRR END----------------------------------------------------
#
#######---------------ChIP_SEQ SRR--------------------------------------------------------
#
#---------------------MACS PEAK CALLING---------------------------------------------------
##This uses MACS to determine ChIP_seq peaks on the individual SRR bam files
#-----------------------------------------------------------------------------------------
MACSpvalue=1e-7
cd ${dir0}/Analysis/ChIP_seq/srr_macs	
##----------------------MACS Peak Calling SRR ChIP seq------------------------------------
for i in $(ls -d */*); do
  cell=$(echo "$i" | cut -d'/' -f1)
  mark=$(echo "$i" | cut -d'/' -f2)
  cd ${dir0}/Analysis/ChIP_seq/srr_macs/$i
  #parallel --xapply --dryrun -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip_srr_${cell}_${mark}.txt "macs14 -t {1} -n {2} -f BAM -g mm -p $MACSpvalue -w --single-profile"
  #parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip_srr_${cell}_${mark}.txt "macs14 -t {1} -n {2} -f BAM -g mm -p $MACSpvalue -w --single-profile"
  #parallel --xapply --dryrun -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip_srr_${cell}_${mark}.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz $datadir/mm10.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig
  #parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip_srr_${cell}_${mark}.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz $datadir/mm10.chrom.sizes {2}.bw" && rm -r *MACS_wiggle #&& rm *negative_peaks*
  #---------------------GENERATE MACS MODEL BUILDING FIGURE---------------------------------
  #for a in `ls *model.r`; do
  #  Rscript --verbose $a
  #done
done
#######---------------ChIP_SEQ SRR END----------------------------------------------------