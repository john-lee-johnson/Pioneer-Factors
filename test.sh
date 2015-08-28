#!/bin/bash
#This script will take ATAC_seq bam files, do MACS peak calling, create bedgraphs of peaks (based on # of tags), and create a file combining all peaks
start=`date +%s`
rm -r *.log
logfile=log.log
exec > $logfile 2>&1

#---------------------Setting Directories-------------------------------------------------
maindir=/mnt/data1/John/Pioneer_Factors
dir0=$1
cd ${dir0}
date > ${dir0}/Time.txt
echo "Start time: $start" >> ${dir0}/Time.txt
#Put all bam files for the cells to be analyzed into a folder named Data

#---------------------MACS PEAK CALLING---------------------------------------------------
##This uses MACS to determine ATAC_seq peaks in the provided bam files
#--------------------OPTIONS--------------------------------------------------------------
option=1	#Set to 1 if you want to run MACS peak calling
#-----------------------------------------------------------------------------------------
MACSpvalue=1e-7
if [ "$option" = "1" ]; then
	end=`date +%s` && runtime=$((end-start)) && echo "Starting MACS Peak Calling" ${runtime} >> ${dir0}/Time.txt

	cd ${dir0}/Analysis/ATAC_seq
	while read line;	##Lists the ATAC seq bam files to do peaking calling
	do
		file="$line"	##ATAC seq bam file
		cell=`echo "$file" | cut -d'_' -f1`		##Cell of interest
		filename="${cell}_ATAC_seq"	##Stores filename
		BAMFILE="${filename}.bam"
		echo "${maindir}/Data/ATAC_seq/bam/${file} ${filename}_macsoutput ${MACSpvalue}"
	done < ${dir0}/options/ATAC_seq_data_files.txt > ${dir0}/options/macs_parallel.txt
	##----------------------MACS Peak Calling--------------------------------------------
	parallel --xapply --dryrun --colsep ' ' -a ${dir0}/options/macs_parallel.txt "macs14 -t {1} -n {2} -p {3} -w --single-profile"
	parallel --xapply --colsep ' ' -a ${dir0}/options/macs_parallel.txt "macs14 -t {1} -n {2} -p {3} -w --single-profile"
	##----------------------Generate BigWig file-----------------------------------------
	parallel --xapply --dryrun --colsep ' ' -a ${dir0}/options/macs_parallel.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz ${dir0}/Data/mm9.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig#parallel --xapply --colsep ' ' -a ${dir0}/options/macs_parallel.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz ${dir0}/Data/mm9.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig
	parallel --xapply --colsep ' ' -a ${dir0}/options/macs_parallel.txt "wigToBigWig -clip {2}_MACS_wiggle/treat/{2}_treat_afterfiting_all.wig.gz ${dir0}/Data/mm9.chrom.sizes {2}_afterfiting_all.bw"	##Converts wig file to BigWig
	rm -r *_macsoutput_MACS_wiggle		##Removes MACS output wig directory
	end=`date +%s` && runtime=$((end-start)) && echo "Time to complete MACS is " ${runtime} >>  ${dir0}/Time.txt
	#Getting total number of tags for ATAC-seq peaks used for normalization later
	cd ${dir0}/Analysis/ATAC_seq
	rm -f ${dir0}/Analysis/ATAC_seq/read_count.txt
	for i in `ls *.xls`;			##Lists MACS output xlx file
	do
		filename=`echo "$i" | cut -d'.' -f1`
		reads=`awk 'BEGIN { FS = ":" } {if (NR == 15) {printf "%.5f\n", $2/1000000}}' ${filename}.xls`			##Gets total number of tags and divides by 10^6 and stores it in an array for later
		cell=`echo $filename | cut -d'_' -f1`
		echo -e "$cell\t$reads"
	done >> ${dir0}/Analysis/ATAC_seq/read_count.txt
	#---------------------CREATE BEDGRAPH FROM PEAKS------------------------------------------
	##This generates a bedgraph of peaks using number of tags from MACS output using the generated xls file
	cd ${dir0}/Analysis/ATAC_seq
	for i in `ls *.xls`;
	do
		filename=`echo "$i" | cut -d'.' -f1`
		${dir0}/Macs2Bedgraph ${filename}.xls > ${filename}.bedgraph	##Custom script to append columns of CHR, Start and End to number of tags from xls file
		sed -i '/_random/d' ${filename}.bedgraph
	done
	#---------------------CREATE TEMP BED FILE------------------------------------------------
	##This generates a temp bed file for unionbedg function later
	cd ${dir0}/Analysis/ATAC_seq
	for i in `ls *peaks.bedgraph`;
	do
		filename=`echo "$i" | cut -d'.' -f1`		 ##Generates filename variable
		tail -n +2 $i > ${filename}_temp.bed 		##Removes the trackline from bedgraph file
	done

	#---------------------UNION FILE----------------------------------------------------------
	#Creating combined peaks file
	#--------------------OPTIONS--------------------------------------------------------------
	cd ${dir0}/Analysis/ATAC_seq
	filename=`ls *_temp.bed  | cut -d'_' -f1 | paste -sd "_" -` 	##Generates combined filename variable
	bedtools unionbedg -i `ls *_temp.bed` | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,5,6,7,8 -o sum,sum,sum,sum,sum > ${filename}_merged_common.bed	##Combines the peak files together and merges bookended peaks and adds the tags for merged peaks to a total sum
	Rscript --verbose ${maindir}/R_Scripts/ATAC_Peaks_Union_Splitter.R $dir0
	rm *_temp.bed		##Removes temporary bed files
fi
#------------------BEDTOOLS COVERAGE------------------------------------------------------
#--------------------OPTIONS--------------------------------------------------------------
option=1	#Set to 1 if you want to run BEDTOOLS COVERAGE
#-----------------------------------------------------------------------------------------
#------------------Sorting Union BED Files to BAM File Order------------------------------
if [ "$option" = "1" ]; then
	end=`date +%s` && runtime=$((end-start)) && echo "Starting bedtools coverage for ATAC seq " ${runtime} >>  ${dir0}/Time.txt
	for i in `ls ${dir0}/Data/ATAC_seq/bam/*.bam`;
	do
		filename=`echo "$i" | rev | cut -d'/' -f1 | rev| cut -d'.' -f1`
		first="`samtools view -H "$i" | grep SQ | cut -f 2 | awk '{ sub(/^SN:/, ""); print;}' | sed -n 1p`"
		if [ "$first" = "chr10" ]; then
			sed -e 's/chr1	/chr1A	/' ${dir0}/Analysis/ATAC_seq/${filename}_union_split.bed | sort -k1,1 -k2,2n | sed -e 's/chr1A/chr1/' > ${dir0}/Analysis/ATAC_seq/${filename}_macsoutput_peaks_sorted.bed
			echo "${dir0}/Analysis/ATAC_seq/${filename}_bam_count.bedgraph ${dir0}/Data/mm9.genome ${dir0}/Analysis/ATAC_seq/${filename}_macsoutput_peaks_sorted.bed ${dir0}/Data/ATAC_seq/bam/${filename}.bam"
		else
			cp ${dir0}/Analysis/ATAC_seq/${filename}_union_split.bed ${dir0}/Analysis/ATAC_seq/${filename}_macsoutput_peaks_sorted.bed
			sed -i '/_random/d' ${dir0}/Analysis/ATAC_seq/${filename}_macsoutput_peaks_sorted.bed
			echo "${dir0}/Analysis/ATAC_seq/${filename}_bam_count.bedgraph ${dir0}/Data/mm9.genome.unsorted ${dir0}/Analysis/ATAC_seq/${filename}_macsoutput_peaks_sorted.bed ${dir0}/Data/ATAC_seq/bam/${filename}.bam"
		fi
	done > ${dir0}/options/ATAC_coverage_parallel.txt
	cd ${dir0}/Analysis/ATAC_seq
	parallel --xapply --dryrun --colsep ' ' -a ${dir0}/options/ATAC_coverage_parallel.txt "bedtools coverage -sorted -g {2} -counts -a {3} -b {4} > {1}"
	#parallel --xapply --colsep ' ' -a ${dir0}/options/ATAC_coverage_parallel.txt "bedtools coverage -sorted -g {2} -counts -a {3} -b {4} > {1}"
	end=`date +%s` && runtime=$((end-start)) && echo "Time to complete coverage of ATAC seq files is " ${runtime} >>  ${dir0}/Time.txt
	cd ${dir0}/Analysis/ATAC_seq
	for i in `ls *_bam_count.bedgraph`;
		do
			filename=`echo "$i" | cut -d'.' -f1`
			awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5}' ${i} > ${filename}_cond.bedgraph	##Removes the old tag counts, keeping only the new tag counts
			sed -e 's/chr1	/chr1A	/' ${filename}_cond.bedgraph | sort -k1,1 -k2,2n | sed -e 's/chr1A/chr1/' > ${filename}_cond_temp.bedgraph
			mv ${filename}_cond_temp.bedgraph ${filename}_cond.bedgraph
		done
	filename="$(ls *bam_count_cond.bedgraph | cut -d'_' -f1 | paste -sd "_" -)"	##Stores a combined filename
	bedtools unionbedg -i `ls *bam_count_cond.bedgraph` | sort -k1,1 -k2,2n | bedtools merge -i stdin -c 4,5,6,7,8 -o sum,sum,sum,sum,sum > ${filename}_bam_count_all.bedgraph
fi
#---------------------R Script that gives ATAC UNIQUE PEAKS and converting to BEDGRAPH----
mkdir -p ${dir0}/Analysis/ATAC_seq/R_output		##Creates R_output directories
mkdir -p ${dir0}/Analysis/ChIP_seq/R_output

Rscript --verbose ${maindir}/R_Scripts/Normalize_ATAC_Peaks.R $dir0
Rscript --verbose ${maindir}/R_Scripts/Normalize_ATAC_Peaks_Percent.R $dir0
Rscript --verbose ${maindir}/R_Scripts/Heatmap1.R $dir0
Rscript --verbose ${maindir}/R_Scripts/Heatmap1_V2.R $dir0
#R --vanilla --slave --args $dir0 < ${maindir}/R_Scripts/Heatmap1.R			##Runs first R script to generate heatmap for ATAC seq unique peaks

#/bin/bash ${dir0}/R_to_Bedgraph.sh			##Converts R output txt files to bedgraph files

#----------------------------ChIP seq analysis--------------------------------------------
#This section of the script will take unique ATAC seq peaks and probe for histone modification data as well as doing motif analysis

#--------------------OPTIONS--------------------------------------------------------------
option=1	#Set to 1 if you want to run bedtools coverage to get histone modification tag counts
#-----------------------------------------------------------------------------------------
if [ "$option" = "1" ]; then
	#------------------Sorting BED Files to BAM File Order------------------------------------
	cd ${dir0}/Analysis/ATAC_seq/R_output
	for i in `ls  *_unique_all.bedgraph`;
	do
		filename=`echo "$i" | cut -d'.' -f1`
		tail -n +2 $i > ${filename}_sorted_temp.bedgraph		##Removes trackline of bedgraph file and stores as temporary file
		sed -e 's/chr1	/chr1A	/' ${filename}_sorted_temp.bedgraph | sort -k1,1 -k2,2n | sed -e 's/chr1A/chr1/' > ${filename}_sorted.bedgraph			##Sorts file
		echo 'track type=bedGraph' | cat - ${filename}_sorted.bedgraph > temp && mv temp ${filename}_sorted.bedgraph		##Adds trackline back
	done
	rm *sorted_temp.bedgraph		#Removes temporary file

	#---------------------Setting Window around each peak-------------------------------------
	##This generates a bedgraph of peaks with a 500bp window around the start and end
	cd ${dir0}/Analysis/ATAC_seq/R_output
	n=1
	for i in `ls  *_unique_all_sorted.bedgraph`;
	do
		name[$((n++))]=`echo "$i" | cut -d'_' -f1`		##Stores names of cells into an array
		filename=`echo "$i" | cut -d'_' -f1`
		tail -n +2 ${dir0}/Analysis/ATAC_seq/R_output/$i | awk 'BEGIN{OFS="\t"} {start=$2-500;end=$3+500;print $1,start,end}' - > ${dir0}/Analysis/ChIP_seq/${filename}_window.bed		##Opens the peaks by 500bp on each side
		mkdir -p ${dir0}/Analysis/ChIP_seq/$filename		##Creates directories for later use
	done
	n=$((n-1))

	#------------COVERAGE HISTONE MODIFICATIONS-----------------------------------------------
	end=`date +%s` && runtime=$((end-start)) && echo "Starting histone mark bedtools coverage " ${runtime} >>  ${dir0}/Time.txt
	while read line; do	##Lists the ATAC seq bam files to do peaking calling
		cell=$line
		for c in `ls ${maindir}/Data/ChIP_seq/{*_B.bam,*_CD4.bam,*_CD8.bam,*_NK.bam,*_HSC.bam,*_MPP.bam,*_CLP.bam}`;
		do
			bamfile=`echo "$c" | rev | cut -d'/' -f1 | rev | cut -d'.' -f1`
			##Gets histone mark tag counts from bam files
			##$cell is the cell of interest and is option -a
			##$c is the bam file to be probed and is option -b
			##Stores in a file named with destination file, then cell (with peaks) listed first, then bam file that was probed
			echo "${dir0}/Analysis/ChIP_seq/$cell/${cell}_${bamfile}_bam_count.bedgraph ${dir0}/Data/mm9.genome ${dir0}/Analysis/ChIP_seq/${cell}_window.bed $c"
		done > ${dir0}/options/${cell}_ChIP_coverage_parallel.txt
		parallel --xapply --dryrun --colsep ' ' -a ${dir0}/options/${cell}_ChIP_coverage_parallel.txt "bedtools coverage -sorted -g {2} -counts -a {3} -b {4} > {1}"
		parallel --xapply --colsep ' ' -a ${dir0}/options/${cell}_ChIP_coverage_parallel.txt "bedtools coverage -sorted -g {2} -counts -a {3} -b {4} > {1}"
		end=`date +%s` && runtime=$((end-start)) && echo "Finished histone mark coverage for ${cell} " ${runtime} >>  ${dir0}/Time.txt
	done < ${dir0}/options/cell.txt
	#Getting total tag counts for all ChIP seq histone modifications
	rm -f ${dir0}/Analysis/ChIP_seq/read_mark.txt		##Removes any old read_mark.txt file
	for i in `ls ${maindir}/Data/ChIP_seq/*.bam`;
	do
		filename=`echo "$i" | cut -d'.' -f1`		##Gets the filename
		mark=`echo "$filename" | cut -d'_' -f1`		##Gets the histone mark
		cell=`echo "$filename" | cut -d'_' -f2,3`			##Gets the cell name
		reads=`samtools idxstats $i | cut -d$'\t' -f3 | paste -sd+ | bc`		##Gets the total number of tags from the bamfile
		#paste <(./progA) <(./progB)
		paste <(echo -e "$mark\t$cell\t$reads") <(echo -e "$reads" | awk '{printf "%.5f\n", $1/1000000}')		##Stores it in a text file
	done >> ${dir0}/Analysis/ChIP_seq/read_mark.txt
	#---------------------UNION---------------------------------------------------------------
	##Creating combined peaks file
	cd ${dir0}/Analysis/ChIP_seq
		for i in `seq 1 $n`;
		do
			cell=${name[i]}		##Gets cell name
			a=1
			while read line; do
			altcell=`echo $line`
			if [ "$cell" != "$altcell" ]; then		##Gets other lymphoid cells
				if [ "$a" == "1" ]; then altcell1=$altcell; fi
				if [ "$a" == "2" ]; then altcell2=$altcell; fi
				if [ "$a" == "3" ]; then altcell3=$altcell; fi
			a=$((a+1))
			fi
		done < ${dir0}/cell.txt
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3}' ${cell}_window.bed > $cell/${cell}_window.bed_temp
	for a in `ls $cell/*count.bedgraph`;
		do
		filename=`echo $a | cut -d'/' -f2 | awk -F'_bam_count.bedgraph' '{print $1}'`
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $4}' $a > $cell/${filename}_temp
		done
		paste `ls $cell/${cell}_window.bed_temp` `ls $cell/${cell}_*_${cell}_temp` `ls $cell/${cell}*LT_HSC_temp` `ls $cell/${cell}*ST_HSC_temp` `ls $cell/${cell}_*_MPP_temp` `ls $cell/${cell}_*_CLP_temp` `ls $cell/${cell}_*_${altcell1}_temp` `ls $cell/${cell}_*_${altcell2}_temp` `ls $cell/${cell}_*_${altcell3}_temp` > ${cell}_merge_bam_count.bedgraph
		echo "CHR" "Start" "End" `ls $cell/${cell}_*_${cell}_temp | cut -d'/' -f2 | cut -d'_' -f2,3` `ls $cell/${cell}_*_LT_HSC_temp | cut -d'/' -f2 | cut -d'_' -f2,3,4` `ls $cell/${cell}_*_ST_HSC_temp | cut -d'/' -f2 | cut -d'_' -f2,3,4` `ls $cell/${cell}_*_MPP_temp | cut -d'/' -f2 | cut -d'_' -f2,3` `ls $cell/${cell}_*_CLP_temp | cut -d'/' -f2 | cut -d'_' -f2,3` `ls $cell/${cell}_*_${altcell1}_temp | cut -d'/' -f2 | cut -d'_' -f2,3` `ls $cell/${cell}_*_${altcell2}_temp | cut -d'/' -f2 | cut -d'_' -f2,3` `ls $cell/${cell}_*_${altcell3}_temp | cut -d'/' -f2 | cut -d'_' -f2,3` | sed -e "s/\ /\t/g" | cat - ${cell}_merge_bam_count.bedgraph > temp && mv temp ${cell}_merge_bam_count.bedgraph               ##Creates a header with column names
		rm $cell/*_temp
	done
	end=`date +%s` && runtime=$((end-start)) && echo "Histone modification coverage files merged, time to complete is " ${runtime} >> ${dir0}/Time.txt

	#---------------NORMALIZATION----------------------------------------------------------
	cd ${dir0}/Analysis/ChIP_seq
	for a in `seq 1 $n`;
	do
		i=`ls ${name[a]}_merge_bam_count.bedgraph`
		filename=`ls ${name[a]}_merge_bam_count.bedgraph | cut -d'.' -f1` 		##Gets filename
		cell=`echo "$filename" | cut -d'_' -f1` 			##Gets name of cell
		b=1
		while read line; do
			altcell=`echo $line`
			if [ "$cell" != "$altcell" ]; then			##Gets other lymphoid cells
				if [ "$b" == "1" ]; then altcell1=$altcell; fi
				if [ "$b" == "2" ]; then altcell2=$altcell; fi
				if [ "$b" == "3" ]; then altcell3=$altcell; fi
				b=$((b+1))
			fi
		done < ${dir0}/cell.txt
		rm -f *temp*
		awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$2,$3}' $i > ${filename}_norm_temp		##Creates a temporary file containing only CHR, Start, End
		while read line; do
			mark="`echo $line | cut -d$' ' -f1`_`echo $line | cut -d' ' -f2`"
			reads=`echo $line | cut -d$' ' -f4`
			awk -v col=$mark -v read=$reads 'BEGIN {FS="\t"} NR==1 {for(i=1;i<=NF;i++){if($i==col){c=i;print $c;break}}} NR>1 {if(i=c){print $c/read}}' $i > ${cell}_${mark}_norm_temp			##Divides the column by the appropriate tag number for normalization and stores as a separate temporary file
		done < read_mark.txt
		find . -size 0 -delete			##Deletes any empty files
		paste `ls ${filename}_norm_temp` `ls ${cell}*_${cell}_norm_temp` `ls ${cell}*_LT_HSC_norm_temp` `ls ${cell}*_ST_HSC_norm_temp` `ls ${cell}*_MPP_norm_temp` `ls ${cell}*_CLP_norm_temp` `ls ${cell}*_${altcell1}_norm_temp` `ls ${cell}*_${altcell2}_norm_temp` `ls ${cell}*_${altcell3}_norm_temp` > ${filename}_norm.bedgraph		##Recombines the temporary files that contain the normalized tag counts into the final file
		cp ${filename}_norm.bedgraph ${filename}_norm_backup.bedgraph
		rm -f *temp*		##Removes the temporary files
	done
fi

#---------------------CREATE 6 COLUMN BED FILE--------------------------------------------
##This generates a 6 column bedgraph that includes unique peak names and adds strand info used for homer
cd ${dir0}/Analysis/ChIP_seq
mkdir -p ${dir0}/Analysis/ChIP_seq/homer
for i in `ls ${dir0}/Analysis/ATAC_seq/R_output/*_unique_all_sorted.bedgraph`;
do
	filename=`echo "$i" | cut -d'.' -f1 | rev | cut -d'/' -f1 | rev`  ##Generates filename variable
	tail -n +2 $i > ${filename}_temp.bed ##Removes the trackline from bedgraph file
	awk 'BEGIN {OFS="\t"} {print $0 "\tpeak" NR "\t0" "\t+"}' ${filename}_temp.bed > ${dir0}/Analysis/ChIP_seq/homer/${filename}_5_column.bed		##Adds a unique peak name and strand info
	rm *_temp.bed 	##Removes the temporary bed file
done

#---------------------R Script that gives histone mark heatmaps---------------------------
R -f ${dir0}/Heatmap2.R		##Runs the second R script that gives heatmaps for histone modifications

/bin/bash /mnt/data1/John/Pioneer-Factors/homer_make_tag.sh
/bin/bash /mnt/data1/John/Pioneer_Factors/peak_annotation.sh
