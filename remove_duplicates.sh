#!/bin/bash

srr=$1
dir=$2

cd $dir
filename=${srr%.*}
rm -f ${filename}_remove.bam
samtools view -f0x400 ${filename} > ${filename}_dups.bam
lines=$(wc -l ${filename}_dups.bam | cut -d' ' -f1)

for ((i=1; i<$lines; i++)); do
line1=$(sed -n "${i}p" < ${filename}_dups.bam)
chr1=$(echo $line1 | cut -d$'\t' -f3 | xargs)
start1=$(echo $line1 | cut -d$'\t' -f4 | xargs)
line2=$(sed -n "${i+1}p" < ${filename}_dups.bam)
chr2=$(echo $line2 | cut -d$'\t' -f3 | xargs)
start2=$(echo $line2 | cut -d$'\t' -f4 | xargs)
line3=$(sed -n "${i+2}p" < ${filename}_dups.bam)
chr3=$(echo $line3 | cut -d$'\t' -f3 | xargs)
start3=$(echo $line3 | cut -d$'\t' -f4 | xargs)
x=1
if [[ "$chr1" = "$chr2" ]] && [[ "$start1" = "$start2" ]]; then
  if [[ "$chr1" = "$chr3" ]] && [[ "$start1" = "$start3" ]]; then
  i=$((i+2))
  linex=$(sed -n "${i+$x}p" < ${filename}_dups.bam)
  chrx=$(echo $linex | cut -d$'\t' -f3 | xargs)
  startx=$(echo $linex | cut -d$'\t' -f4 | xargs)
    while [[ "$chr1" = "$chrx" ]] && [[ "$start1" = "$startx" ]]; do
    echo "$linex" >> ${filename}_remove.bam
    echo "successfully removed line"
    x=$((x+1))
    linex=$(sed -n "${i+$x}p" < ${filename}_dups.bam)
    chrx=$(echo $linex | cut -d$'\t' -f3 | xargs)
    startx=$(echo $linex | cut -d$'\t' -f4 | xargs)
    done
    i=$(($i + $x))
  fi
fi
done