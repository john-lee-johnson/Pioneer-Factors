#!/bin/bash
#Produces normalized bigwig files


##---------------------NORMALIZE BIG WIG FILES FROM BAM-----------------------------------
#Produces normalized bigwig files from the aligned bam file
lines=$(wc -l $infodir/sample_files_amit.txt | cut -d' ' -f1)
for ((i=1; i<$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_amit.txt)
gsm=$(echo $line | cut -d':' -f1 | xargs)
seq=$(echo $line | cut -d':' -f2 | cut -d';' -f3 | xargs)
cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq/$cell
  filename=${cell}_${seq}_Amit_mm10_${gsm}
  fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/Amit/$filename/tagInfo.txt" | cut -d"=" -f2)
  echo "Generating bigWig for $filename"
  bamCoverage --bam ${filename}.bam --outFileName ${filename}.bw --outFileFormat bigwig --bamIndex ${filename}.bam.bai --normalizeUsingRPKM --binSize 10 --centerReads --binSize 10 --numberOfProcessors 38 --fragmentLength $fragLength
  cp ${filename}.bw $datadir/ATAC_seq/bw
fi
if [[ "$seq" = ChIP* ]]; then
  mark=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs | cut -d'_' -f1 | xargs)
  cell=$(echo $line | cut -d':' -f2 | cut -d';' -f1 | xargs | cut -d'_' -f2 | xargs)
  filename=${cell}_${seq}_${mark}_Amit_mm10_${gsm}
  cd $dir0/Data/ChIP_seq/$cell/$mark
  fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/Amit/$filename/tagInfo.txt" | cut -d"=" -f2)
  echo "Generating bigWig for $filename"
  bamCoverage --bam ${filename}.bam --outFileName ${filename}.bw --outFileFormat bigwig --bamIndex ${filename}.bam.bai --normalizeUsingRPKM --binSize 10 --centerReads --binSize 10 --numberOfProcessors 38 --fragmentLength $fragLength
  cp ${filename}.bw $datadir/ChIP_seq/bw
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq/$cell
  filename=${cell}_${seq}_Amit_mm10_${gsm}
fi
done