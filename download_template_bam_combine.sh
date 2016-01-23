#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory

##---------------------SETTING FUNCTIONS FOR COMBINING REPLICATES-------------------------
#Combine bam files for technical replicates
function bam_combine {
unset input
unset reads
for a in $(ls SRR*.bam | grep -v STAR | grep -v adapter | grep -v marked | grep -v removed | grep -v 'sort' | grep -v reads); do
  input=$(echo $input" INPUT="`pwd`/$a)
done
echo "PICARD VERSION MergeSamFiles: `java -jar $PICARD MergeSamFiles --version`" >> $logdir/download_log_${des}.txt
echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MergeSamFiles $input ""OUTPUT="`pwd`/${filename}.bam >> $paralleldir/merge_bam_${des}.txt
}

rm -f $paralleldir/merge_bam_${des}.txt

echo""; echo "COMBINING BAM FILES FOR TECHNICAL REPLICATES"
echo""; echo "#--------------COMBINING BAM FILES--------------------"
lines=$(wc -l $infodir/sample_files_${des}.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}.txt)
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
investigator=$(echo $line | cut -d',' -f5)
sequencer=$(echo $line | cut -d',' -f6)
type=$(echo $line | cut -d',' -f7)
replicate=$(echo $line | cut -d',' -f8)
rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell/$investigator ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  echo "Combining Bam files for $filename"
  bam_combine
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  investigator=$(echo $line | cut -d',' -f6)
  sequencer=$(echo $line | cut -d',' -f7)
  type=$(echo $line | cut -d',' -f8)
  replicate=$(echo $line | cut -d',' -f9)
  rep=$(echo $line | cut -d',' -f9 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f9 | cut -d'_' -f2)
  cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
  echo "Combining Bam files for $filename"
  bam_combine
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  echo "Combining Bam files for $filename"
  bam_combine
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  echo "Combining Bam files for $filename"
  bam_combine
fi
done
echo "#--------------------------------------------------------"
echo ""; echo "To Combine Bam Files, set combine_bam=true and combine_fastq=true"; echo ""

echo "--------COMBINING BAM FILES--------"; echo ""
if [ -f $paralleldir/merge_bam_${des}.txt ]; then
  echo "#----------------COMBINING TECHNICAL REPLICATES-----------------------" >> $logdir/download_log_${des}.txt
  parallel --xapply --dryrun -j 30 -- < $paralleldir/merge_bam_${des}.txt >> $logdir/download_log_${des}.txt
  parallel --xapply --dryrun -j 30 -- < $paralleldir/merge_bam_${des}.txt
  parallel --xapply -j 30 -- < $paralleldir/merge_bam_${des}.txt
fi