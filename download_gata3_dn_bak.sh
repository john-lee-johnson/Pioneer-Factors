#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject
#Process and align the files
#Combine technical and biological replicates
#Do peak calling
#Move files to the data directory

#If the SAMPLE FILE is properly FORMATTED, set test1=true
test1=true
#If the SRR FILE is properly FORMATTED, set test2=true
test2=true
#To CARRY OUT the download, set download=true
#IF already DOWNLOADED, set download=false
download=false
#If you have CHECKED for quality and adapter contamination and set the options correctly, set test3=true to CARRY OUT preprocessing
#Set options for preprocessing below
test3=true
#To CARRY OUT the preprocessing, set process=true
#If files already PROCESSED, set process=false
process=false
#To ALIGN, set align = true
#If ALIGNMENT has already been carried out, set align=false
align=false
#To remove DUPLICATES, set duplicate=true
#If remove DUPLICATES has already been carried out, set duplicate=false
duplicate=false
     #Set removed=true to remove duplicates
     #Set removed=false to mark duplicates, but keep them
     removed=true
#When files have been processed, aligned, deduplicated, and are now ready to be combined, set test4=true
test4=true
#To CARRY OUT the actual combining of bam files, set combine_bam=true
#IF files already COMBINED, set combine_bam=false
combine_bam=false
#To CARRY OUT the actual combining of fastq files, set combine_fastq=true
#IF files already COMBINED, set combine_fastq=false
combine_fastq=false
#To CARRY OUT the combining of biological replicates, set replicate=true
#IF files already COMBINED, set replicate=false
replicate=false
#To CARRYOUT Homer makeTagDirectory, set homer=true
#IF tag directories already created, set homer=false
homer=false
#For MACS14 Peak Calling, set macs_14=true
#IF Peak Calling is already done, set macs_14=false
macs_14=false
#Check to make sure the MACS14 commands file is correct before setting an option and carrying out the peak calling
macs_14_chip=false #For TF ChIP-seq (with Input) MACS14 set macs_14_chip=true
macs_14_atac=false #For ATAC-seq MACS14 set macs_14_chip=true
macs_14_model=false #To build Model PDF from R script, set macs_14_model=true
macs_14_bigwig=false #To normalize BigWig file from MACS, set macs_14_bigwig=true
macs_14_bedgraph=false #To convert MACS14 XLS to Bedgraph and do Bedtools Coverage, set macs_14_bedgraph=true
#To copy files to data dir, set copy = true
#If files already COPIED, set copy = false
copy_bam=false


##---------------------SET BIOPROJECT NUMBER----------------------------------------------
query=PRJNA124395

##---------------------SET GSE NUMBER HERE------------------------------------------------
#GSE number will identify the info and parallel command files being output in this script
des=GSE20898

#--------------------PREPROCESSING OPTIONS------------------------------------------------
#Set commands for quality trimming, adapter trimming here:
#Set the proper adapter
atacquality=false
atacadapter=false
adapteratac="nextera"
chipquality=false
chipadapter=false
adapterchip="illumina"
rnaquality=false
rnaadapter=false
adapterrna="small_rna"
mnasequality=false
mnaseadapter=false
#--illumina              Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter
#--nextera               Adapter sequence to be trimmed is the first 12bp of the Nextera adapter
#--small_rna             Adapter sequence to be trimmed is the first 12bp of the Illumina Small RNA Adapter
#-----------------------------------------------------------------------------------------

#Sets the working directories
infodir=/mnt/data1/John/Pioneer-Factors/info
maindir=/mnt/data1/John/Pioneer_Factors
scriptdir=/mnt/data1/John/Pioneer-Factors
filedir=/mnt/data1/John/Pioneer-Factors/sample_files
paralleldir=/mnt/data1/John/Pioneer-Factors/parallel_commands
r_dir=/mnt/data1/John/Pioneer-Factors/R_Scripts
datadir=/mnt/data1/VahediLab/PTF_Team/Data
datadir2=/mnt/data1/VahediLab/PTF_Team/Data2
dir0=$maindir

rm -f $infodir/srr_files_${des}.txt
rm -f $paralleldir/srr_download_${des}.txt
rm -f $paralleldir/process_fastq_${des}.txt
rm -f $paralleldir/remove_duplicates_${des}.txt
rm -f $paralleldir/mark_duplicates_${des}.txt
rm -f $paralleldir/merge_bam_${des}.txt
rm -f $paralleldir/homer_key_file.txt
rm -f $paralleldir/macs_parallel_${des}.txt

##---------------------SETTING FUNCTIONS FOR DOWNLOAD-------------------------------------
#Download functions should account for proper header information and paired vs. single end data

#Creates parallel commands for fastq-dump single-end Illumina data
function fastq-download {
  fastq-dump -F -X 1 -Z ${srr} > $infodir/test.txt #Will download one line of the file to check for header information
    if [[ $(head $infodir/test.txt | sed -n 1p | cut -d":" -f1) = "@1" ]]; then
      echo "fastq-dump ${srr} -O ${wd}" >> $paralleldir/srr_download_${des}.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "fastq-dump -F ${srr} -O ${wd}"  >> $paralleldir/srr_download_${des}.txt #If file does have proper header information, will keep the original header
    fi
}
#Creates parallel commands for fastq-dump paired-end Illumina data
function fastq-download-paired {
  fastq-dump -F -X 1 -Z ${srr} > $infodir/test.txt #Will download one line of the file to check for header information
    if [[ $(head $infodir/test.txt | sed -n 1p | cut -d":" -f1) = "@1" ]]; then
      echo "fastq-dump -I --split-files ${srr} -O ${wd}" >> $paralleldir/srr_download_${des}.txt #If file does not have proper header information, will keep the SRR header
    else
      echo "fastq-dump -I --split-files -F ${srr} -O ${wd}"  >> $paralleldir/srr_download_${des}.txt #If file does have proper header information, will keep the original header
    fi
}
#Creates parallel commands for fastq-dump single-end ABI data
function abi-download {
    echo "abi-dump -F ${srr} -O ${wd}" >> $paralleldir/srr_download_${des}.txt #If file does not have proper header information, will keep the SRR header
}
#Creates parallel commands for fastq-dump paired-end ABI data
function abi-download-paired {
    echo "abi-dump -F ${srr} -O ${wd}" >> $paralleldir/srr_download_${des}.txt #If file does not have proper header information, will keep the SRR header
}

##---------------------SETTING FUNCTIONS FOR ALIGNMENT-------------------------------------
#Alignment functions should ultimately generate the aligned file in the following format:
#SRR123456789.bam

#STAR generates genome index
function star_generate {
STAR --runMode genomeGenerate --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --genomeFastaFiles /mnt/data0/John/GRCm38p4_mm10_pa_only.fa --sjdbGTFfile /mnt/data0/John/gencode.vM6.annotation.gtf
}
#STAR alignment for ATAC-seq
function star_atac {
STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate 
mv Aligned.sortedByCoord.out.bam STAR_${srr}.bam #Keeps the bam file that STAR generates
mv Log.final.out ${srr}_Log.final.out #Keeps a copy of the aligner log
STAR --runMode inputAlignmentsFromBAM --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --inputBAMfile STAR_${srr}.bam --bamRemoveDuplicatesType UniqueIdentical
mv Processed.out.bam dedupSTAR_${srr}.bam #Keeps a bam file with PCR duplicates removed via STAR
}
#ABI alignment using BOWTIE
function abi-align {
cp /mnt/data0/John/bowtie_mm10_colorspace/mm10_colorspace.*.ebwt $wd
mkdir -p tmp
bowtie -p 38 -S -C mm10_colorspace -f ${srr}_F3.csfasta -Q ${srr}_F3_QV.qual ${srr}.sam
samtools view -b -h ${srr}.sam > ${srr}_all_reads.bam
rm -f ${srr}.sam
samtools view -b -h -F 4 ${srr}_all_reads.bam > ${srr}.bam
mv ${srr}.bam ${srr}_sort.bam
echo "------------SORTING BAM FILE-----------------"
java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD SortSam INPUT=${srr}_sort.bam OUTPUT=${srr}.bam SORT_ORDER=coordinate
#cp ${srr}_F3.csfasta $datadir/$seq/reads
#cp ${srr}_F3_QV.qual $datadir/$seq/reads
}

#STAR alignment for ChIP-seq
function star_chip {
STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate 
mv Aligned.sortedByCoord.out.bam ${srr}.bam #Keeps the bam file that STAR generates
mv Log.final.out ${srr}_Log.final.out #Keeps a copy of the aligner log
}

#STAR alignment for RNA-seq 2-pass
function star_rna {
STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30
mv Aligned.sortedByCoord.out.bam ${srr}.bam #Keeps the bam file that STAR generates
mv Log.final.out ${srr}_Log.final.out #Keeps a copy of the aligner log
STAR --runMode inputAlignmentsFromBAM --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --inputBAMfile ${srr}.bam --bamRemoveDuplicatesType UniqueIdentical
mv Processed.out.bam dedupSTAR_${srr}.bam #Keeps a bam file with PCR duplicates removed via STAR
}

##---------------------SETTING FUNCTIONS FOR DUPLICATE REMOVAL----------------------------
#Duplicate removal should account for proper header information
#A backup file should be kept of duplicates that are marked, but not removed

function mark_duplicates {
if [[ $(head -1 ${srr}.fastq | cut -c 1-4) = "@SRR" ]]; then
  #If file does not have header information regarding location of the read on the lane, will not look for optical duplicates
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_removed.bam READ_NAME_REGEX=null REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_${des}.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_marked.bam READ_NAME_REGEX=null METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_${des}.txt
else
  #If file has proper header information regarding location of the read on the lane, will look for optical duplicates
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_removed.bam REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_${des}.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_marked.bam METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_${des}.txt
fi
}

function mark_duplicates_abi {
  #Will not look for optical duplicates because ABI SOLiD 
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_removed.bam READ_NAME_REGEX=null REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_${des}.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_marked.bam READ_NAME_REGEX=null METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_${des}.txt
  #java -Xmx2g -jar $PICARD MarkDuplicates INPUT=`pwd`/${srr}.bam OUTPUT=`pwd`/${srr}_removed.bam READ_NAME_REGEX=null REMOVE_DUPLICATES=true METRICS_FILE=`pwd`/${srr}_metrics.txt
  #java -Xmx2g -jar $PICARD MarkDuplicates INPUT=`pwd`/${srr}.bam OUTPUT=`pwd`/${srr}_marked.bam READ_NAME_REGEX=null METRICS_FILE=`pwd`/dupsMarked_${srr}_metrics.txt
}

#STAR remove genome from memory
function star_remove {
STAR --genomeLoad Remove --genomeDir /mnt/data0/John/genome_GRCm38p4_M6
}

#Processes fastq files to remove adapters
function adapter-trim {
  #echo "trim_galore -o "`pwd`" --${adapter} "`pwd`"/${srr}.fastq" >> $paralleldir/process_fastq_${des}.txt
  cp ${srr}.fastq ${srr}_adapter.fastq
  trim_galore -o `pwd` --${adapter} `pwd`/${srr}.fastq
  echo "Adapter trimming: trim_galore -o `pwd` --${adapter} `pwd`/${srr}.fastq"
  mv ${srr}_trimmed.fq ${srr}.fastq
}

#Processes fastq files to remove low quality reads
function quality-trim {
  echo "Quality Trimming"
  #echo "fastq_masker -q 20 -i `pwd`/${srr}_trimmed.fq -o `pwd`/${srr}_mask.fastq" >> $paralleldir/process_fastq_${des}.txt #Will mask low quality data
  cp ${srr}.fastq ${srr}_quality.fastq
  fastq_masker -q 20 -i `pwd`/${srr}.fastq -o `pwd`/${srr}_mask.fastq #Will mask low quality data
  mv ${srr}_mask.fastq ${srr}.fastq
}

#Combine bam files for replicates
function bam_combine {
unset input
unset reads
for a in $(ls SRR*.bam | grep -v STAR | grep -v adapter | grep -v marked | grep -v removed | grep -v 'sort' | grep -v reads); do
  input=$(echo $input" INPUT="`pwd`/$a)
done
echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MergeSamFiles $input ""OUTPUT="`pwd`/${filename}.bam >> $paralleldir/merge_bam_${des}.txt
}

#Combine fastq files for replicates
function fastq_combine {
allFastq=$(ls SRR*.fastq | grep -v STAR | grep -v adapter | grep -v marked | grep -v removed | grep -v 'sort' | grep -v reads)
cat $allFastq > ${filename}.fastq
pigz -f -p 35 ${filename}.fastq
}

function replicate-combine {
unset input
unset reads
for a in $(ls *GSM*.bam | grep -v STAR | grep -v adapter | grep -v marked | grep -v removed | grep -v 'sort' | grep -v reads); do
  input=$(echo $input" INPUT="`pwd`/$a)
done
echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MergeSamFiles $input ""OUTPUT="`pwd`/${filename}.bam >> $paralleldir/merge_bam_${des}.txt
allFastq=$(ls *GSM*.fastq.gz | grep -v STAR | grep -v adapter | grep -v marked | grep -v removed | grep -v 'sort' | grep -v reads)
cat $allFastq > ${filename}.fastq.gz
}

##---------------------DOWNLOAD BIOPROJECT SAMPLE DATA------------------------------------
esearch -db sra -query $query | efetch --format docsum | xtract -pattern DocumentSummary -element Title > $infodir/sample_files_${des}.txt


##---------------------Notes about Format-------------------------------------------------
#The format for the sample file needs to be in the following format:
#gsm,cell,species,seq,(mark),sequencer,type,replicate
#gsm should look like GSM101236
#cell should look like DP
#species should look like 'Mus musculus'
#seq should look like 'ChIP_seq'
#If seq is histone modification ChIP-seq, seq should look like ChIP_seq,H3K4me1
#If seq is transcription factor ChIP-seq, seq should look like ChIP_seq,TCF1
#If seq is input DNA ChIP-seq, seq should look like ChIP_seq,Input
#sequencer should be either ABI or Illumina
#type should be either Paired or Single
#If there are biological replicates, replicate should look like Rep_1 or Rep_2
#If there are not biological replicates, replicate should look liks Rep_0

##---------------------CLEAN UP SAMPLE FILE-----------------------------------------------
#Here are sample sed commands:
#Replace text with NewText
###sed -i -e 's/text/NewText/g' $infodir/sample_files_${des}.txt
#Delete lines with pattern
###sed -i -e '/pattern/d' $infodir/sample_files_${des}.txt
#Print lines matching pattern
###sed -i -n '/pattern;/p' $infodir/sample_files_${des}.txt
#sed -i -e 's/text/NewText/g' $infodir/sample_files_${des}.txt
#sed -i -e 's/text/NewText/g' $infodir/sample_files_${des}.txt
#sed -i -e 's/text/NewText/g' $infodir/sample_files_${des}.txt
#sed -i -n '/pattern;/p' $infodir/sample_files_${des}.txt
#sed -i -n '/ThyDP/p' $infodir/sample_files_${des}.txt
sed -i -e '/Methyl-Seq/d' $infodir/sample_files_${des}.txt
sed -i -e '/RNA-seq/d' $infodir/sample_files_${des}.txt
sed -i -e '/Th2-/d' $infodir/sample_files_${des}.txt
sed -i -e '/CD8-/d' $infodir/sample_files_${des}.txt
sed -i -e '/NKT-/d' $infodir/sample_files_${des}.txt
sed -i -e '/nTreg/d' $infodir/sample_files_${des}.txt
sed -i -e '/iTreg/d' $infodir/sample_files_${des}.txt
sed -i -e '/Th17/d' $infodir/sample_files_${des}.txt
sed -i -e '/Th1/d' $infodir/sample_files_${des}.txt
sed -i -e '/CD4/d' $infodir/sample_files_${des}.txt
sed -i -e 's/\[ChIP-seq\]/ChIP-seq/g' $infodir/sample_files_${des}.txt
sed -i -e 's/\[ChIP-Seq\]/ChIP-Seq/g' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM523222: DP-Gata3 ChIP-seq/GSM523222,DP,Mus musculus,ChIP_seq,Gata3_Zhao,Illumina,Single,Rep_1/g' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM523221: DN-Gata3 ChIP-seq/GSM523221,DN,Mus musculus,ChIP_seq,Gata3_Zhao,Illumina,Single,Rep_0/g' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM742021: DP-Gata3-replicate ChIP-Seq/GSM742021,DP,Mus musculus,ChIP_seq,Gata3_Zhao,Illumina,Single,Rep_2/g' $infodir/sample_files_${des}.txt

echo "#--------------------SAMPLE FILES-----------------------"
sed '' $infodir/sample_files_${des}.txt
echo "#--------------------------------------------------------"
echo "If the SAMPLE FILE is properly formated, set test1 = true"

if [[ $test1 = true ]]; then
echo ""; echo "----MOVING ONTO DOWNLOAD PREP-----"

##---------------------DOWNLOAD SRR NUMBERS ----------------------------------------------
#Will read in GSM numbers from earlier and download the SRR numbers
#SRR File should have the following format:
#Example_seq (Mark) Cell Dir SRR1 SRR2 SRR3 ... 
echo ""
echo "#--------------------SRR FILES---------------------------"
lines=$(wc -l $infodir/sample_files_${des}.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}.txt)
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; mkdir -p $cell ; cd $cell ; if [[ "$repnum" -ne 0 ]]; then mkdir -p $replicate; cd $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq at `pwd`"
  echo "$gsm,$cell,$species,$seq,$sequencer,$type,$replicate,`pwd`:"`esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  cd $dir0/Data/ChIP_seq ; mkdir -p $cell/$mark ; cd $cell/$mark ; if [[ "$repnum" -ne 0 ]]; then mkdir -p $replicate; cd $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq $mark at `pwd`"
  echo "$gsm,$cell,$species,$seq,$mark,$sequencer,$type,$replicate,`pwd`:"`esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; mkdir -p $cell ; cd $cell ; if [[ "$repnum" -ne 0 ]]; then mkdir -p $replicate; cd $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq at `pwd`"
  echo "$gsm,$cell,$species,$seq,$sequencer,$type,$replicate,`pwd`:"`esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; mkdir -p $cell ; cd $cell ; if [[ "$repnum" -ne 0 ]]; then mkdir -p $replicate; cd $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq at `pwd`"
  echo "$gsm,$cell,$species,$seq,$sequencer,$type,$replicate,`pwd`:"`esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
done
echo "#--------------------------------------------------------"
echo "SRR File should have the following format gsm,cell,species,seq,(mark),sequencer,type,pwd:SRR1 SRR2 SRR3 ..." 
echo ""
sed '' $infodir/srr_files_${des}.txt
echo "If the SRR FILE is properly formated, set test2 = true"; echo ""

if [[ $test2 = true ]]; then
echo "----MOVING ONTO ACTUAL DOWNLOAD-----"

###---------------------SRR DATA PARALLEL COMMANDS OUTPUT---------------------------------
#Read in a file and set up parallel commands for later SRR download
echo "--------------------FASTQ-DUMP--------------------"
while read line; do
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
wd=$(echo $line | cut -d',' -f8 | cut -d':' -f1)
if [[ "$seq" = ATAC* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "download Paired Illumina data"
        fastq-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "download Single-end Illumina data"
        fastq-download
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "download Paired ABI data"
        abi-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "download Single-end ABI data"
        abi-download
      fi
    fi
  done
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  wd=$(echo $line | cut -d',' -f9 | cut -d':' -f1)
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "download"
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        fastq-download-paired
      elif [[ "$type" = "Single" ]]; then
        fastq-download
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        abi-download-paired
      elif [[ "$type" = "Single" ]]; then
        abi-download
      fi
    fi
  done
fi
if [[ "$seq" = RNA* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "download"
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        fastq-download-paired
      elif [[ "$type" = "Single" ]]; then
        fastq-download
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        abi-download-paired
      elif [[ "$type" = "Single" ]]; then
        abi-download
      fi
    fi
  done
fi
if [[ "$seq" = MNase* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "download"
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        fastq-download-paired
      elif [[ "$type" = "Single" ]]; then
        fastq-download
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        abi-download-paired
      elif [[ "$type" = "Single" ]]; then
        abi-download
      fi
    fi
  done
fi
done < $infodir/srr_files_${des}.txt
echo "#--------------------------------------------------------"
echo ""
echo "-----------------CHECK PARALLEL COMMANDS-----------------"
sed '' $paralleldir/srr_download_${des}.txt
echo "#--------------------------------------------------------"
echo ""

echo "Check to make sure download commands are correct before setting download = true"
if [[ $download = true ]]; then
parallel --xapply --dryrun -j 30 -- < $paralleldir/srr_download_${des}.txt
parallel --xapply -j 30 -- < $paralleldir/srr_download_${des}.txt
fi #download

echo ""; echo "Recommended to stop here and check a few of the files for quality and adapter contamination"
echo "Set test3 = true when ready to continue"; echo ""
                        
if [[ $test3 = true ]]; then
echo "MOVING ONTO PREPROCESSING"; echo ""
echo "--------------------Process Files--------------------"
if [[ $process = true ]]; then
echo "PROCESSING BAM FILES"; echo ""
while read line; do
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
wd=$(echo $line | cut -d',' -f8 | cut -d':' -f1)
if [[ "$seq" = ATAC* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "quality or adapter trim"
    if [[ $atacquality = true ]]; then
      quality-trim
    fi
    if [[ $atacadapter = true ]]; then
      adapter=$adapteratac
      adapter-trim
    fi
  done
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  wd=$(echo $line | cut -d',' -f9 | cut -d':' -f1)
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "quality or adapter trim"
    if [[ $chipadapter = true ]]; then
      adapter=$adapterchip
      adapter-trim
    fi
    if [[ $chipquality = true ]]; then
      echo "Skipping quality trimming"
      #quality-trim
    fi
  done
fi
if [[ "$seq" = RNA* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "quality or adapter trim"
    if [[ $rnaquality = true ]]; then
      quality-trim
    fi
    if [[ $rnaadapter = true ]]; then
      adapter=$adapterrna
      adapter-trim
    fi
  done
fi
if [[ "$seq" = MNase* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "quality or adapter trim"
    if [[ $mnasequality = true ]]; then
      quality-trim
    fi
    if [[ $mnaseadapter = true ]]; then
      adapter=$adaptermnase
      adapter-trim
    fi
  done
fi
done < $infodir/srr_files_${des}.txt
fi #process

echo ""; echo "Once files are downloaded and preprocessed, set align = TRUE to align SRR files"
if [[ $align = true ]]; then
###---------------------ALIGNMENT---------------------------------
#Read in SRR files and align
echo "--------------------ALIGNMENT--------------------"
while read line; do
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
wd=$(echo $line | cut -d',' -f8 | cut -d':' -f1)
if [[ "$seq" = ATAC* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "align"
    #SET ALIGNMENT FUNCTION HERE----------------------------
  done
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  wd=$(echo $line | cut -d',' -f9 | cut -d':' -f1)
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "align"
    star_chip
  done
fi
if [[ "$seq" = RNA* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "align"
    #SET ALIGNMENT FUNCTION HERE----------------------------
  done
fi
if [[ "$seq" = MNase* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "align"
    #SET ALIGNMENT FUNCTION HERE----------------------------
  done
fi
done < $infodir/srr_files_${des}.txt
echo "#--------------------------------------------------------"
echo ""
fi #align

echo ""; echo "Once files have been aligned, set duplicate = TRUE to mark and remove duplicates"; echo ""
if [[ $duplicate = true ]]; then
echo "---------------REMOVE DUPLICATES--------------------"
while read line; do
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
wd=$(echo $line | cut -d',' -f8 | cut -d':' -f1)
if [[ "$seq" = ATAC* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "duplicate removal"
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates_abi
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates_abi
      fi
    fi
  done
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  wd=$(echo $line | cut -d',' -f9 | cut -d':' -f1)
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "duplicate removal"
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates_abi
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates_abi
      fi
    fi
  done
fi
if [[ "$seq" = RNA* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "duplicate removal"
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates_abi
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates_abi
      fi
    fi
  done
fi
if [[ "$seq" = MNase* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "duplicate removal"
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        mark_duplicates_abi
      elif [[ "$type" = "Single" ]]; then
        mark_duplicates_abi
      fi
    fi
  done
fi
done < $infodir/srr_files_${des}.txt

if [[ $removed = true ]]; then
  if [ -f $paralleldir/remove_duplicates_${des}.txt ]; then
  echo "Removing duplicates"
  parallel --xapply --dryrun -j 30 -- < $paralleldir/remove_duplicates_${des}.txt
  parallel --xapply -j 30 -- < $paralleldir/remove_duplicates_${des}.txt
fi
else
if [ -f $paralleldir/mark_duplicates_${des}.txt ]; then
  echo "Marking, but not removing duplicates"
  parallel --xapply --dryrun -j 30 -- < $paralleldir/mark_duplicates_${des}.txt
  parallel --xapply -j 30 -- < $paralleldir/mark_duplicates_${des}.txt
fi
fi

echo""; echo "---------------COPY DEDUPLICATED FILES--------------------"
while read line; do
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
wd=$(echo $line | cut -d',' -f8 | cut -d':' -f1)
if [[ "$seq" = ATAC* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ $removed = true ]]; then
      echo $srr "Removing Duplicated Sequences"
      cp ${srr}_removed.bam ${srr}.bam
    else
      echo $srr "Keeping Duplicated Sequences"
      cp ${srr}_marked.bam ${srr}.bam
    fi
  done
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  wd=$(echo $line | cut -d',' -f9 | cut -d':' -f1)
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ $removed = true ]]; then
      echo $srr "Removing Duplicated Sequences"
      cp ${srr}_removed.bam ${srr}.bam
    else
      echo $srr "Keeping Duplicated Sequences"
      cp ${srr}_marked.bam ${srr}.bam
    fi
  done
fi
if [[ "$seq" = RNA* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ $removed = true ]]; then
      echo $srr "Removing Duplicated Sequences"
      cp ${srr}_removed.bam ${srr}.bam
    else
      echo $srr "Keeping Duplicated Sequences"
      cp ${srr}_marked.bam ${srr}.bam
    fi
  done
fi
if [[ "$seq" = MNase* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f8 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ $removed = true ]]; then
      echo $srr "Removing Duplicated Sequences"
      cp ${srr}_removed.bam ${srr}.bam
    else
      echo $srr "Keeping Duplicated Sequences"
      cp ${srr}_marked.bam ${srr}.bam
    fi
  done
fi
done < $infodir/srr_files_${des}.txt
fi #duplicate

if [[ $test4 = true ]]; then
echo""; echo "MOVING ONTO COMBINING BAM FILES"
echo""; echo "---------------COMBINING BAM FILES--------------------"
lines=$(wc -l $infodir/sample_files_${des}.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}.txt)
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_mm10_${gsm}
  echo "Combining Bam files for $filename"
  bam_combine
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  cd $dir0/Data/ChIP_seq ; cd $cell/$mark ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_${mark}_mm10_${gsm}
  echo "Combining Bam files for $filename"
  bam_combine
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_mm10_${gsm}
  echo "Combining Bam files for $filename"
  bam_combine
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_mm10_${gsm}
  echo "Combining Bam files for $filename"
  bam_combine
fi
done
echo "#--------------------------------------------------------"
echo ""; echo "To Combine Bam Files, set combine_bam=TRUE"; echo ""

if [[ $combine_bam = true ]]; then
echo "COMBINING BAM FILES"; echo ""
if [ -f $paralleldir/merge_bam_${des}.txt ]; then
parallel --xapply --dryrun -j 30 -- < $paralleldir/merge_bam_${des}.txt
parallel --xapply -j 30 -- < $paralleldir/merge_bam_${des}.txt
fi
fi #combine

if [[ $combine_fastq = true ]]; then
echo""; echo "MOVING ONTO COMBINING FASTQ FILES"
echo""; echo "---------------COMBINING FASTQ FILES--------------------"
lines=$(wc -l $infodir/sample_files_${des}.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}.txt)
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_mm10_${gsm}
  echo "Combining fastq files for $filename"
  fastq_combine
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  cd $dir0/Data/ChIP_seq ; cd $cell/$mark ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_${mark}_mm10_${gsm}
  echo "Combining fastq files for $filename"
  fastq_combine
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_mm10_${gsm}
  echo "Combining fastq files for $filename"
  fastq_combine
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_mm10_${gsm}
  echo "Combining fastq files for $filename"
  fastq_combine
fi
done
echo "#--------------------------------------------------------"
echo ""; echo "To Combine Files, set combine = TRUE"; echo ""
fi #combine_fastq

if [[ $replicate = true ]]; then
echo "Combining Biological Replicates"; echo ""
echo "#--------------------Combining Biological Replicates---------------------------"
echo "#-------------------Part 1 Copying Files---------------------------"
lines=$(wc -l $infodir/sample_files_${des}.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}.txt)
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell
  if [[ "$repnum" -ne 0 ]]; then 
    cd $replicate
    filename=${cell}_${seq}_mm10_${gsm}
    echo "Copying Rep $repnum to upper directory"
    cp ${filename}.bam ..
    cp ${filename}.fastq ..
  fi
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  cd $dir0/Data/ChIP_seq ; cd $cell/$mark
  if [[ "$repnum" -ne 0 ]]; then
    cd $replicate
    filename=${cell}_${seq}_${mark}_mm10_${gsm}
    echo "Copying Rep $repnum to upper directory"
    cp ${filename}.bam ..
    cp ${filename}.fastq ..
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell
  if [[ "$repnum" -ne 0 ]]; then
    cd $replicate
    filename=${cell}_${seq}_mm10_${gsm}
    echo "Copying Rep $repnum to upper directory"
    cp ${filename}.bam ..
    cp ${filename}.fastq ..
  fi
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell
  if [[ "$repnum" -ne 0 ]]; then
    cd $replicate
    filename=${cell}_${seq}_mm10_${gsm}
    echo "Copying Rep $repnum to upper directory"
    cp ${filename}.bam ..
    cp ${filename}.fastq ..
  fi
fi
done
fi #replicate
rm -f $paralleldir/merge_bam_${des}.txt
rm -f $infodir/sample_files_${des}_noreps.txt

echo "#--------------Making New Sample File and Homer Key File---------------------------"
lines=$(wc -l $infodir/sample_files_${des}.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}.txt)
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell
  if [[ "$repnum" = 0 ]]; then 
    filename=${cell}_${seq}_mm10_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_mm10_${gsm}
    replicate-combine
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  cd $dir0/Data/ChIP_seq ; cd $cell/$mark
  if [[ "$repnum" = 0 ]]; then 
    filename=${cell}_${seq}_${mark}_mm10_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${mark}_mm10_${gsm}
    replicate-combine
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell
  if [[ "$repnum" = 0 ]]; then 
  filename=${cell}_${seq}_mm10_${gsm}
  echo "$line" >> $infodir/sample_files_${des}_noreps.txt
  echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_mm10_${gsm}
    replicate-combine
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell
  if [[ "$repnum" = 0 ]]; then 
  filename=${cell}_${seq}_mm10_${gsm}
  echo "$line" >> $infodir/sample_files_${des}_noreps.txt
  echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_mm10_${gsm}
    replicate-combine
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
fi
done
echo "#--------------------------------------------------------"
echo ""
sed '' $infodir/sample_files_${des}_noreps.txt

if [[ $replicate = true ]]; then
echo "COMBINING BAM FILES"; echo ""
if [ -f $paralleldir/merge_bam_${des}.txt ]; then
parallel --xapply --dryrun -j 30 -- < $paralleldir/merge_bam_${des}.txt
parallel --xapply -j 30 -- < $paralleldir/merge_bam_${des}.txt
fi
fi #replicate

if [[ $homer = true ]]; then
echo "MAKING HOMER TAG DIRECTORIES"; echo ""
if [ -f $paralleldir/homer_key_file.txt ]; then
cd $dir0/Analysis/Homer/Tag_Directories/; mkdir -p $des ; cd $des
sed '' $paralleldir/homer_key_file.txt
batchMakeTagDirectory.pl $paralleldir/homer_key_file.txt -cpu 40 -genome mm10 -format sam
fi
fi #homer

if [[ $macs_14 = true ]]; then
##--------------MACS 14 OUTPUT COMMANDS COMBINED BAM FILE-------------------------------
#Read in a file and output commands to a file for MACS peak calling 14
echo "#-------------------Part 1 Copying Files---------------------------"
lines=$(wc -l $infodir/sample_files_${des}_noreps.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}_noreps.txt)
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell ; mkdir -p macs
  filename=${cell}_${seq}_mm10_${gsm}
  fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/$des/$filename/tagInfo.txt" | cut -d"=" -f2 | xargs)
  echo "`pwd`/${filename}.bam `pwd`/macs/${filename}" "$fragLength">> $paralleldir/macs_parallel_${des}.txt
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  filename=${cell}_${seq}_${mark}_mm10_${gsm}
  if echo "$mark" | grep -q Input; then
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark
    echo "Input found: `pwd`/${filename}.bam"
  else
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark ; mkdir -p macs
    fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/$des/$filename/tagInfo.txt" | cut -d"=" -f2 | xargs)
    echo "`pwd`/${filename}.bam `pwd`/macs/${filename}" "$fragLength">> $paralleldir/macs_parallel_${des}.txt
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
done

##---------------------CLEAN UP MACS14 PARALLEL COMMANDS-----------------------------------------------
#Here are sample sed commands:
#Replace text with NewText
###sed -i -e 's/text/NewText/g' $paralleldir/macs_parallel_${des}.txt
#Delete lines with pattern
###sed -i -e '/pattern/d' $paralleldir/macs_parallel_${des}.txt
#Print lines matching pattern
###sed -i -n '/pattern;/p' $paralleldir/macs_parallel_${des}.txt
sed -i -e 's, 65, 65 /mnt/data1/John/Pioneer_Factors/Data/ChIP_seq/DP/Input_Rothenberg/DP_ChIP_seq_Input_Rothenberg_mm10_GSM774302.bam,g' $paralleldir/macs_parallel_${des}.txt
sed -i -e 's, 210, 210 /mnt/data1/John/Pioneer_Factors/Data/ChIP_seq/DP/Input_Rothenberg/DP_ChIP_seq_Input_Rothenberg_mm10_GSM774302.bam,g' $paralleldir/macs_parallel_${des}.txt
#sed -i -n '/pattern;/p' $paralleldir/macs_parallel_${des}.txt

echo "#--------------------------------------------------------"
echo ""
echo "-----------------CHECK MACS14 PARALLEL COMMANDS-----------------"
sed '' $paralleldir/macs_parallel_${des}.txt
echo "#--------------------------------------------------------"
echo ""
echo "Check to make sure MACS14 commands are correct before setting the next option"
echo "Format should be PATH_TO_BAM OUTPUT_NAME fragLength Input"
echo "Transcription factor ChIP-seq should have the full path to the Input file added to the end of each line (separated with a space)"
echo "Set macs_14_chip = true if sample files are transcription factor ChIP-seq with an Input file"
echo "Set macs_14_atac = true if sample files do not have an Input file"

#---------Clean up MACS Parallel commands-------------------------

#---------------------MACS 14 PEAK CALLING---------------------------------------------------
##This uses MACS to determine ATAC_seq peaks on the replicate combined bam files
#-----------------------------------------------------------------------------------------
MACSpvalue=1e-7

##----------------------MACS Peak Calling ChIP seq----------------------------------------
if [[ $macs_14_chip = true ]]; then
echo ""; echo "--------------MACS Peak Calling ChIP seq---------------"
if [ -f $paralleldir/macs_parallel_${des}.txt ]; then
parallel --dryrun --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_${des}.txt "macs14 -t {1} -n {2} -c {4} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"
parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_${des}.txt "macs14 -t {1} -n {2} -c {4} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"
fi
fi #macs_14_chip 

##----------------------MACS Peak Calling ATAC seq----------------------------------------
if [[ $macs_14_atac = true ]]; then
echo ""; echo "--------------MACS Peak Calling ATAC seq---------------"
if [ -f $paralleldir/macs_parallel_${des}.txt ]; then
parallel --xapply --dryrun -j 35 --colsep ' ' -a $paralleldir/macs_parallel_${des}.txt "macs14 -t {1} -n {2} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"
parallel --xapply -j 35 --colsep ' ' -a $paralleldir/macs_parallel_${des}.txt "macs14 -t {1} -n {2} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"
fi
fi #macs_14_atac

if [[ $macs_14_model = true ]]; then
#---------------------GENERATE MACS MODEL BUILDING FIGURE---------------------------------
#Generate PDF of Model Building
lines=$(wc -l $infodir/sample_files_${des}_noreps.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}_noreps.txt)
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell ; cd macs
  filename=${cell}_${seq}_mm10_${gsm}
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  filename=${cell}_${seq}_${mark}_mm10_${gsm}
  if echo "$mark" | grep -q Input; then
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark
    echo "Skipping Input"
  else
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark ; cd macs
    for a in `ls *model.r`; do
      echo "Generating MACS model figure PDF for ${filename}"
      Rscript --verbose $a
    done
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell
  echo "No Peak Calling for RNA-seq"
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell
  echo "No Peak Calling for MNase-seq"
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
done
fi #macs_14_model

#---------------------CHECK MACS ERRORS---------------------------------------------------
#Not in this current version

#---------------------Normalize WIG FILE--------------------------------------------------
echo ""; echo "-------------Normalize WIG FILE---------------"
if [[ $macs_14_bigwig = true ]]; then
lines=$(wc -l $infodir/sample_files_${des}_noreps.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}_noreps.txt)
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell ; cd macs
  filename=${cell}_${seq}_mm10_${gsm}
  wd=`pwd`
  if [ ! -f ${filename}_normalized.wig ]; then
  if [ -f ${filename}_treat_afterfiting_all.wig.gz ]; then
    echo "$filename normalizing of wig file underway"
    gzip -d ${filename}_treat_afterfiting_all.wig.gz
    readsTotal=$(samtools view -c ../${filename}.bam)
    reads=$(bc <<< "scale=2;$readsTotal/1000000")
    echo "$reads"
    perl $scriptdir/spmr.pl "$reads" ${filename}_treat_afterfiting_all.wig ${filename}_normalized.wig
    wigToBigWig ${filename}_normalized.wig $datadir/mm10.chrom.sizes ${filename}_normalized.bw -clip
    rm -r ${filename}_MACS_wiggle
  fi
  else
     echo "$filename Normalization already completed"
  fi
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  filename=${cell}_${seq}_${mark}_mm10_${gsm}
  if echo "$mark" | grep -q Input; then
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark
    echo "Skipping Input"
  else
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark ; cd macs
      wd=`pwd`
      if [ ! -f ${filename}_normalized.wig ]; then
      if [ -f ${filename}_treat_afterfiting_all.wig.gz ]; then
        echo "$filename normalizing of wig file underway"
        gzip -d ${filename}_treat_afterfiting_all.wig.gz
        readsTotal=$(samtools view -c ../${filename}.bam)
        reads=$(bc <<< "scale=2;$readsTotal/1000000")
        echo "$reads"
        perl $scriptdir/spmr.pl "$reads" ${filename}_treat_afterfiting_all.wig ${filename}_normalized.wig
        wigToBigWig ${filename}_normalized.wig $datadir/mm10.chrom.sizes ${filename}_normalized.bw -clip
        rm -r ${filename}_MACS_wiggle
        rm -f ${filename}_control_afterfiting_all.wig.gz
      fi
      else
        echo "$filename Normalization already completed"
      fi
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell
  echo "No Peak Calling for RNA-seq"
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell
  echo "No Peak Calling for MNase-seq"
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
done
fi #macs_14_bigwig

#------------------Make Bedgraph from Bedtools Coverage-----------------------------------
if [[ $macs_14_bedgraph = true ]]; then
echo ""; echo "-------------------Make Bedgraph from Bedtools Coverage----------------------"
lines=$(wc -l $infodir/sample_files_${des}_noreps.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}_noreps.txt)
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell ; cd macs
  filename=${cell}_${seq}_mm10_${gsm}
  wd=`pwd`
  #---------------------CREATE BEDGRAPH FROM PEAKS------------------------------------------
  ##This generates a bedgraph of peaks using number of tags from MACS output using the generated xls file 
  $scriptdir/Macs2Bedgraph ${filename}_peaks.xls > ${filename}.bedgraph	##Custom script to append columns of CHR, Start and End to number of tags from xls file
  #sed -i '/_random/d' ${filename}.bedgraph
  sed -e 's/chr1	/chr1A	/' ${filename}.bedgraph  | sort -k1,1 -k2,2n | sed -e 's/chr1A/chr1/' > ${filename}_sort.bedgraph 
  bedtools coverage -sorted -g $datadir/sorted.mm10.genome -counts -a ${filename}_sort.bedgraph  -b ../${filename}.bam | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' - > ${filename}_tag_peaks.bedgraph
  rm -f ${filename}.bedgraph
  rm -f ${filename}_sort.bedgraph
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  filename=${cell}_${seq}_${mark}_mm10_${gsm}
  if echo "$mark" | grep -q Input; then
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark
    echo "Skipping Input"
  else
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark ; cd macs
    echo "Generating bedgraph for $filename"
    wd=`pwd`
    #---------------------CREATE BEDGRAPH FROM PEAKS------------------------------------------
    ##This generates a bedgraph of peaks using number of tags from MACS output using the generated xls file 
    $scriptdir/Macs2Bedgraph ${filename}_peaks.xls > ${filename}.bedgraph	##Custom script to append columns of CHR, Start and End to number of tags from xls file
    #sed -i '/_random/d' ${filename}.bedgraph
    sed -e 's/chr1	/chr1A	/' ${filename}.bedgraph  | sort -k1,1 -k2,2n | sed -e 's/chr1A/chr1/' > ${filename}_sort.bedgraph 
    bedtools coverage -sorted -g $datadir/sorted.mm10.genome -counts -a ${filename}_sort.bedgraph  -b ../${filename}.bam | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' - > ${filename}_tag_peaks.bedgraph
    rm -f ${filename}.bedgraph
    rm -f ${filename}_sort.bedgraph
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell
  echo "No Peak Calling for RNA-seq"
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell
  echo "No Peak Calling for MNase-seq"
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
done
fi #macs_14_bedgraph

fi #macs_14

if [[ $copy_bam = true ]]; then
echo "COPYING BAM FILES TO DATA DIR"; echo ""
echo "#--------------------Copy Bam File to Data Dir---------------------------"
lines=$(wc -l $infodir/sample_files_${des}_noreps.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}_noreps.txt)
gsm=$(echo $line | cut -d',' -f1)
cell=$(echo $line | cut -d',' -f2)
species=$(echo $line | cut -d',' -f3)
seq=$(echo $line | cut -d',' -f4)
sequencer=$(echo $line | cut -d',' -f5)
type=$(echo $line | cut -d',' -f6)
replicate=$(echo $line | cut -d',' -f7)
rep=$(echo $line | cut -d',' -f7 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f7 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell
  filename=${cell}_${seq}_mm10_${gsm}
  echo "Copying Bam and MACS files for $filename to Data Directory"
  samtools index ${filename}.bam
  cp ${filename}.bam $datadir2/ATAC_seq/bam
  cp ${filename}.bam.bai $datadir2/ATAC_seq/bam
  cp -R macs $datadir2/ATAC_seq/
  cp macs/${filename}_normalized.bw $datadir2/ChIP_seq/bw
  cp ${filename}.fastq.gz $datadir2/ChIP_seq/fastq
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  sequencer=$(echo $line | cut -d',' -f6)
  type=$(echo $line | cut -d',' -f7)
  replicate=$(echo $line | cut -d',' -f8)
  rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
  cd $dir0/Data/ChIP_seq ; cd $cell/$mark
  filename=${cell}_${seq}_${mark}_mm10_${gsm}
  if echo "$mark" | grep -q Input; then
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark
    echo "Skipping Input"
  else
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark
    echo "Copying Bam and MACS files for $filename to Data Directory"
    samtools index ${filename}.bam
    cp ${filename}.bam $datadir2/ChIP_seq/bam
    cp ${filename}.bam.bai $datadir2/ChIP_seq/bam
    cp -R macs $datadir2/ChIP_seq/
    cp macs/${filename}_normalized.bw $datadir2/ChIP_seq/bw
    cp ${filename}.fastq.gz $datadir2/ChIP_seq/fastq
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell
  filename=${cell}_${seq}_mm10_${gsm}
  echo "Copying Bam files for $filename to Data Directory"
  samtools index ${filename}.bam
  cp ${filename}.bam $datadir2/RNA_seq/bam
  cp ${filename}.bam.bai $datadir2/RNA_seq/bam
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell
  filename=${cell}_${seq}_mm10_${gsm}
  echo "Copying Bam files for $filename to Data Directory"
  cp ${filename}.bam $datadir2/MNase_seq/bam
  cp ${filename}.bam.bai $datadir2/MNase_seq/bam
fi
done
echo "#--------------------------------------------------------"
echo ""; echo "Finished copying"; echo ""
fi #copy_bam

fi #test4
fi #test3
fi #test2
fi #test1
