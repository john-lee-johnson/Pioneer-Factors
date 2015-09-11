#!/bin/bash
for i in `ls -d */`;
do
	filename=`echo "$i" | rev | cut -d'/' -f2 | rev`		 ##Generates filename variable
	cd $filename
	fastq=`ls *.fastq`
#	STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/mm9.M6.index --readFilesIn $fastq --outSAMtype BAM SortedByCoordinate
	mv Aligned.sortedByCoord.out.bam $filename.bam
	cd ..
	#tail -n +2 $i > ${filename}_temp.bed 		##Removes the trackline from bedgraph file
done
