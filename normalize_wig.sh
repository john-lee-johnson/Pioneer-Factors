#!/bin/bash
cd $wd

#---------------------Normalize WIG FILE--------------------------------------------------
for i in `ls -d *MACS_wiggle*`; do
#  echo $i
  cd $i/treat
  filename=${i%_macs_MACS_wiggle}
  if [ ! -f ${filename}_macs_normalized.wig ]; then
  if [ -f ${filename}_macs_treat_afterfiting_all.wig.gz ]; then
    echo "$filename normalizing of wig file underway"
    gzip -d ${filename}_macs_treat_afterfiting_all.wig.gz
     bamDir=$(find $datadir -name "${filename}.bam")
     readsTotal=$(samtools view -c $bamDir)
     reads=$(bc <<< "scale=2;$readsTotal/1000000")
     perl $scriptdir/spmr.pl "$reads" ${filename}_macs_treat_afterfiting_all.wig ${filename}_macs_normalized.wig
  fi
  else
     echo "$filename Normalization already completed"
  fi
  cd $wd
done
