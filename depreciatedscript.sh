#####-------BROAD PEAK CALLING DEPRECIATED------------------------------------------------
cd ${dir0}/Analysis/ChIP_seq/macs
for i in `ls {DN3.*peaks.bed,DP.*peaks.bed}`; do
filename=`echo "$i" | rev | cut -d'.' -f2,3 | rev`
awk 'BEGIN{ OFS="\t"}{printf "%s\t%s\t%s\t%s\t%.5f\n", $1, $2, $3, $4, $5/10}' $i > ${filename}_test.bed
done


####----------------REPLICATES------------------------------------------------------------
mkdir -p $dir0/Analysis/ChIP_seq/mspc
cd $dir0/Analysis/ChIP_seq/mspc
#rm -r Session_01
#mono /mnt/data1/bin/MSPCCLI_v1_Mono/MSPC.exe -i ${dir0}/Analysis/ChIP_seq/macs/DN3.1_H3Ac_macs_14_peaks_test.bed -i ${dir0}/Analysis/ChIP_seq/macs/DN3.2_H3Ac_macs_14_peaks_test.bed -r biological -s 1E-7 -w 1E-4 -g 1E-8
cp Session_01/*DN3.1_H3Ac_macs_14_peaks_test/A_Output_Set.bed ./DN3.1_Output_Set.bed
cp Session_01/*DN3.2_H3Ac_macs_14_peaks_test/A_Output_Set.bed ./DN3.2_Output_Set.bed
for i in `ls  *Output_Set.bed`; do
  filename=`echo "$i" | rev | cut -d'_' -f2,3 | rev`
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5}' $i > ${filename}.bedgraph
done
bedtools unionbedg -i DN3.1_Output.bedgraph DN3.2_Output.bedgraph > Output_Set_Union.bed 
bedtools merge -i Output_Set_Union.bed > Output_Set_Merge.bedgraph


#Getting total number of tags for ATAC-seq peaks used for normalization later
#cd ${dir0}/Data/ChIP_seq/tcf1
#for i in `ls Thy/Tcf1/Log.final.out`; do
#cell=$(echo $i | cut -d'/' -f1)
#echo -e $cell $(awk 'BEGIN { FS = "|\t" } {if (NR == 9) {printf "%.5f\n", $2/1000000}}' $i)
#done
##----------------------------------------------------------------------------------------
##--------------------------------END TCF1 END--------------------------------------------

	
#	#Getting total number of tags for ATAC-seq peaks used for normalization later
#	cd ${dir0}/Analysis/ATAC_seq
#	rm -f ${dir0}/Analysis/ATAC_seq/read_count.txt
#	for i in `ls *.xls`;			##Lists MACS output xlx file
#	do
#		filename=`echo "$i" | cut -d'.' -f1`
#		reads=`awk 'BEGIN { FS = ":" } {if (NR == 15) {printf "%.5f\n", $2/1000000}}' ${filename}.xls`			##Gets total number of tags and divides by 10^6 and stores it in an array for later
#		cell=`echo $filename | cut -d'_' -f1`
#		echo -e "$cell\t$reads"
#	done >> ${dir0}/Analysis/ATAC_seq/read_count.txt
