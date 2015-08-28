#!/bin/bash
#---------------------Setting Directories-------------------------------------------------
maindir=/mnt/data1/John/Pioneer_Factors
dir0=/mnt/data0/John/Pioneer_Factors/homer
cd ${dir0}
#---------------------HOMER PEAK ANNOTATION-----------------------------------------------
cd ${maindir}/Analysis/ChIP_seq/homer
mkdir -p ${dir0}/peak_annotation
mkdir -p ${dir0}/heatmap
for i in `ls *5_column.bed`;
do
	maincell=`echo "$i" | cut -d'_' -f1`			##The cell of interest
	mkdir -p ${dir0}/peak_annotation/$maincell
	mkdir -p ${dir0}/heatmap/$maincell
	for c in `ls -d ${dir0}/ChIP_seq/H3*/`;
	do
		mark=`echo ${c} | rev | cut -d'/' -f2 | rev | cut -d'_' -f1`	##Gets the histone marks
		cell=`echo ${c} | rev | cut -d'/' -f2 | rev | cut -d'_' -f2,3`	##The cell to be probed
		echo -e "$i $c $maincell $cell $mark"
		#annotatePeaks.pl ${i} mm9 -size 5000 -hist 50 -d ${c} > ${maincell}_$cell_${c2}_${c1}_output.txt	##Gets quantitative info for histone marks at $maincell unique peaks
		#annotatePeaks.pl ${i} mm9 -size 5000 -hist 50 -ghist -d ${c} > ${maincell}_$cell_${c2}_${c1}_heatmap.txt		##Gets heatmap for histone marks at $maincell unique peaks
	done
done > ${dir0}/peak_annotation.txt
parallel --xapply --dryrun --colsep ' ' -a ${dir0}/peak_annotation.txt "annotatePeaks.pl {1} mm9 -size 5000 -hist 50 -d {2} > ${dir0}/peak_annotation/{3]/{3}_{4}_{5}_output.txt"
parallel --xapply --colsep ' ' -a ${dir0}/peak_annotation.txt "annotatePeaks.pl {1} mm9 -size 5000 -hist 50 -d {2} > ${dir0}/peak_annotation/{3}/{3}_{4}_{5}_output.txt"
parallel --xapply --dryrun --colsep ' ' -a ${dir0}/peak_annotation.txt "annotatePeaks.pl {1} mm9 -size 5000 -hist 50 -ghist -d {2} > ${dir0}/heatmap/{3]/{3}_{4}_{5}_heatmap.txt"
parallel --xapply --colsep ' ' -a ${dir0}/peak_annotation.txt "annotatePeaks.pl {1} mm9 -size 5000 -hist 50 -ghist -d {2} > ${dir0}/heatmap/{3}/{3}_{4}_{5}_heatmap.txt"
