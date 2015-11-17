#!/bin/bash
rm -f $paralleldir/macs_parallel_redo.txt
rm -f $paralleldir/macs_parallel_redo_bigwig.txt
cd $wd	
for i in `ls *.log`; do
filename=${i%_macs.log}
line=$(sed -n "8p" < $i)
line2=$(sed -n '/predicted fragment length/p' $i)
line3=$(sed -n -e 's/^.*predicted fragment length is //p' <$i)
echo "$filename $line compared to $line3"
done

for i in `ls *.log`; do
  if grep -Fq WARNING $i;
  then
    filename=${i%_macs.log}
    fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/Amit/$filename/tagInfo.txt" | cut -d"=" -f2 | xargs)
    echo "$filename $fragLength"
    bamDir=$(find $datadir -name "${filename}.bam")
    echo "$bamDir"
    echo "macs14 -t $bamDir -n ${filename}_macs --bw $fragLength -f BAM -g mm -p 1e-7 -m 5,30 -w --single-profile >> $filename_macs.log 2>&1" >> $paralleldir/macs_parallel_redo.txt
    echo "wigToBigWig -clip ${filename}_macs_MACS_wiggle/treat/${filename}_macs_treat_afterfiting_all.wig.gz /mnt/data1/VahediLab/PTF_Team/Data/mm10.chrom.sizes ${filename}_macs_afterfiting_all.bw" >> $paralleldir/macs_parallel_redo_bigwig.txt
    mv ${filename}_macs.log ${filename}_macs.log.old
    mv ${filename}_macs.pdf ${filename}_macs_old.pdf
    rm -r ${filename}*MACS_wiggle*
  fi
done

if [ -f $paralleldir/macs_parallel_redo.txt ]; then
parallel --xapply --dryrun -j 40 -- < $paralleldir/macs_parallel_redo.txt
parallel --xapply -j 40 -- < $paralleldir/macs_parallel_redo.txt
fi
if [ -f $paralleldir/macs_parallel_redo_bigwig.txt ]; then
parallel --xapply --dryrun -j 40 -- < $paralleldir/macs_parallel_redo_bigwig.txt
parallel --xapply -j 40 -- < $paralleldir/macs_parallel_redo_bigwig.txt
fi

for i in `ls *model.r`; do
  Rscript --verbose $i
done
