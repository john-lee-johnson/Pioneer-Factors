#!/bin/bash
#---------------------Setting Directories-------------------------------------------------
maindir=/mnt/data1/John/Pioneer_Factors
dir0=/mnt/data0/John/Pioneer_Factors/homer
cd ${dir0}
#---------------------HOMER PEAK ANNOTATION-----------------------------------------------
cd ${maindir}/Analysis/ChIP_seq/homer
for i in `ls *5_column.bed`;
do
	maincell=`echo "$i" | cut -d'_' -f1`			##The cell of interest
	for n in `ls *5_column.bed`;
	do
		cell=`echo "$n" | cut -d'_' -f1`			##The cell of interest
			for c in `ls -d ${dir0}/ChIP_seq/H3*/`;
			do
				mark=`echo ${c} | rev | cut -d'/' -f3 | rev | cut -d'_' -f1`		##Gets the histone marks
				echo -e "$i\t$n\t$mark"
				#annotatePeaks.pl ${i} mm9 -size 5000 -hist 50 -d ${c} > ${maincell}_$cell_${c2}_${c1}_output.txt	##Gets quantitative info for histone marks at $maincell unique peaks
				#annotatePeaks.pl ${i} mm9 -size 5000 -hist 50 -ghist -d ${c} > ${maincell}_$cell_${c2}_${c1}_heatmap.txt		##Gets heatmap for histone marks at $maincell unique peaks
			done
	done
done
