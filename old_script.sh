


##---------------------DOWNLOAD SRR AND ALIGN TF and OUTPUT COMMANDS----------------------
#Read in a file and download the SRR file and convert to fastq for TF Data
rm -f $paralleldir/macs_parallel_tf.txt
while read line; do
  name=$(echo $line | cut -d' ' -f1)
  gsm=$(echo $line | cut -d' ' -f2)
  srr=$(echo $line | cut -d' ' -f3)
  cell=$(echo $line | cut -d' ' -f4 | cut -d'_' -f1)
  tf=$(echo $line | cut -d' ' -f4 | cut -d'_' -f2)
  filename1=${cell}_${tf}_${name}_mm10_${gsm}
  fastq=${filename1}.fastq
  cd $dir0/Data/ChIP_seq/Transcription_Factor ; mkdir -p $tf ; cd $tf
  echo $filename1
  #-------------FASTQ DUMP and STAR TF----------------------------------------------------
  #fastq-dump $srr && cp ${srr}.fastq ${filename1}.fastq
  #pigz -f -p 40  ${filename1}.fastq
  #star_chip 
  #samtools view -F0x400 -b Processed.out.bam > ${filename1}.bam
  #fastqc ${filename1}.fastq.gz
  #echo $tf $cell $(samtools view -c ${filename1}.bam) $(samtools view -c Processed.out.bam) >> $infodir/pcr_duplicates.txt
  #&& mv Aligned.sortedByCoord.out.bam ${filename1}.bam
  #mv ${filename1}.bam $datadir/ChIP_seq/bam
  #samtools index $datadir/ChIP_seq/bam/${filename1}.bam
  #mv ${filename1}.fastq.gz $datadir/ChIP_seq/fastq
  gsm=$(echo $line | cut -d' ' -f5)
  srr=$(echo $line | cut -d' ' -f6)
  cell=$(echo $line | cut -d' ' -f7 | cut -d'_' -f1)
  tf=$(echo $line | cut -d' ' -f7 | cut -d'_' -f2)
  filename2=${cell}_${tf}_${name}_mm10_${gsm}
  fastq=${filename2}.fastq
  mkdir -p Input ; cd Input
  echo $filename2
  #-------------FASTQ DUMP and STAR INPUT-------------------------------------------------
  #fastq-dump $srr && cp ${srr}.fastq ${filename2}.fastq #&& rm ${srr}.fastq
  #pigz -f -p 40  ${filename2}.fastq
  #star_chip #&& mv Aligned.sortedByCoord.out.bam ${filename2}.bam
  #samtools view -F0x400 -b Processed.out.bam > ${filename2}.bam
  #fastqc ${filename2}.fastq.gz
  #echo $tf $cell $(samtools view -c ${filename2}.bam) $(samtools view -c Processed.out.bam) >> $infodir/pcr_duplicates.txt
  #mv ${filename2}.bam $datadir/ChIP_seq/bam
  #samtools index $datadir/ChIP_seq/bam/${filename2}.bam
  #mv ${filename2}.fastq.gz $datadir/ChIP_seq/fastq
  #---------------------MACS OUTPUT COMMANDS TF--------------------------------------------
  echo "$datadir/ChIP_seq/bam/${filename1}.bam ${filename1}_macs $datadir/ChIP_seq/bam/${filename2}.bam" >> $paralleldir/macs_parallel_tf.txt
done < $filedir/sample_transcription_factors.txt

##---------------------DOWNLOAD SRR ROTHENBERG--------------------------------------------
#Read in a file and download the SRR file and convert to fastq for Amit Data
while read line; do
  srr=$(echo $line | cut -d' ' -f1)
  seq=$(echo $line | cut -d' ' -f2 | cut -d'_' -f1)
  cell=$(echo $line | cut -d' ' -f2 | cut -d'_' -f2)
if [[ "$seq" = H3* ]]
then
  echo $srr $seq $cell
  cd $dir0/Data/ChIP_seq ; mkdir -p $cell ; cd $cell ; mkdir -p $seq
  cd $seq
  #fastq-download
fi
if [[ "$seq" = RNA ]]
then
  echo $srr $seq $cell
  cd $dir0/Data/RNA_seq ; mkdir -p $cell ; cd $cell
  #fastq-dump $srr
fi
done < $filedir/srr_files_rothenberg.txt

##---------------------COMBINE AND ALIGN ROTHENBERG---------------------------------------
#Combine fastq files, compress to .gz, and align
cd $dir0
while read line; do
  gsm=$(echo $line | cut -d' ' -f1)
  seq=$(echo $line | cut -d' ' -f2)
  cell=$(echo $line | cut -d' ' -f3)
  filename=${cell}_${seq}_Rothenberg_mm10_${gsm}
  fastq=${filename}.fastq
  echo $seq $cell
if [[ "$seq" = H3* ]]; then
  cd $dir0/Data/ChIP_seq/$cell/$seq
  #cat SRR*.fastq > ${filename}.fastq
  #pigz -p 40  ${fastq} #&& rm SRR*.fastq
  #star_chip 
  #samtools view -F0x400 -b Processed.out.bam > ${filename}.bam
  #fastqc ${filename}.fastq.gz
  #echo $seq $cell $(samtools view -c ${filename}.bam) $(samtools view -c Processed.out.bam) >> $infodir/pcr_duplicates.txt
  #&& mv Aligned.sortedByCoord.out.bam ${filename}.bam
  #mv ${filename}.bam $datadir/ChIP_seq/bam
  #samtools index $datadir/ChIP_seq/bam/${filename}.bam
  #mv ${fastq}.gz $datadir/ChIP_seq/fastq
fi
cd $dir0
if [[ "$seq" = RNA ]]; then
  cd $dir0/Data/RNA_seq/$cell
  #cat SRR*.fastq > ${filename}.fastq
  #pigz -p 40  ${fastq} #&& rm SRR*.fastq
  #star_rna && mv Aligned.sortedByCoord.out.bam ${filename}.bam
  #mv ${filename}.bam $datadir/RNA_seq/bam
  #samtools index $datadir/RNA_seq/bam/${filename}.bam
  #mv ${fastq}.gz $datadir/RNA_seq/fastq
fi
done < $filedir/sample_files_rothenberg.txt

##---------------------MACS OUTPUT COMMANDS ROTHENBERG------------------------------------
#Read in a file and output commands to a file for MACS peak calling for Amit Data
cd $dir0
rm -f $paralleldir/macs_parallel_chip_rothenberg.txt
while read line; do
  gsm=$(echo $line | cut -d' ' -f1)
  seq=$(echo $line | cut -d' ' -f2)
  cell=$(echo $line | cut -d' ' -f3)
  filename=${cell}_${seq}_Rothenberg_mm10_${gsm}
  fastq=${filename}.fastq
cd $dir0
if [[ "$seq" = H3* ]] 
then
  cd $dir0/Data/ChIP_seq/$cell/$seq
  echo "$datadir/ChIP_seq/bam/${filename}.bam ${filename}_macs" >> $paralleldir/macs_parallel_chip_rothenberg.txt
fi
cd $dir0
if [[ "$seq" = RNA ]]
then
  cd $dir0/Data/RNA_seq/$cell
fi
cd $dir0
done < $filedir/sample_files_rothenberg.txt

#---Remove STAR Genome from Memory--------------------------------------------------------
cd $dir0
star_remove
rm -r _STARtmp
rm Log.out
rm Log.progress.out
rm Aligned.out.sam

cd ${paralleldir}
split -l 39 macs_parallel_chip.txt macs_parallel_chip
