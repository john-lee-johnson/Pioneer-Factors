#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject
#Process and align the files
#Combine technical and biological replicates
#Do peak calling
#Generate normalized bigwig files
#Move files to the data directory

#-------------------------HOW TO and OPTIONS----------------------------------------------
#Step 1: Add the BIOPROJECT NUMBER of the sample to be downloaded
#Step 2: Add the GSE NUMBER of the sample to be downloaded 
#Step 3: Run the script, use SED to properly format the file according to guidelines
#If the SAMPLE FILE is properly FORMATTED, set test1=true
test1=true
#Step 4: Check to make sure the SRR FILE is properly formatted
#If the SRR FILE is properly FORMATTED, set test2=true
test2=true
#To CARRY OUT the download, set download=true
#Once DOWNLOADED, set download=false
download=true
#Step 5: Use FASTQC to check the fastq files for quality and adapter contaminations
#Step 6: Set the options for preprocessing (below) and then set test3=true to continue
#Set options for preprocessing below
test3=true
#To CARRY OUT the preprocessing, set process=true
#Once files are PROCESSED, set process=false
process=true
#To ALIGN, set align = true
#Once ALIGNED, set align=false
align=true
#To remove DUPLICATES, set duplicate=true
#Once DUPLICATES have been removed, set duplicate=false
#Note: If you want to keep duplicates, set removed=false (below) before running duplicate=true
duplicate=true
#Step 7: When files have been processed, aligned, deduplicated, and technical replicates are now ready to be combined, set test4=true
test4=true
#To CARRY OUT the actual combining of bam files, set combine_bam=true
#Once files are COMBINED, set combine_bam=false
combine_bam=true
#To CARRY OUT the actual combining of fastq files, set combine_fastq=true
#IF files already COMBINED, set combine_fastq=false
combine_fastq=true
#Step 8: To combine biological replicates, set replicate_copy=true to copy the files for the replicates into the right place
#Once files are put into the right place, set replicate_copy=false
replicate_copy=true
#Once bam and fastq files combined for technical replicates are combined, and biological replicates have been copied, set test5=true
test5=true
#To COMBINE biological replicates, set replicate_combine=true
#Once COMBINED, set replicate_combine=false
replicate_combine=true
#To CARRYOUT Homer makeTagDirectory, set homer=true
#Once tag directories are CREATED, set homer=false
homer=true
#Step 9: Choose or omit Peak Calling options below
#For MACS14 Peak Calling, set macs_14=true
macs_14=true
#Check to make sure the MACS14 commands file is correct (below) before setting an option and CARRYING OUT the peak calling
macs_14_chip=true #For TF ChIP-seq (with Input) MACS14 set macs_14_chip=true
macs_14_atac=false #For ATAC-seq MACS14 set macs_14_chip=true
macs_14_model=true #To build Model PDF from R script, set macs_14_model=true
macs_14_bigwig=true #To normalize BigWig file from MACS, set macs_14_bigwig=true
macs_14_bedgraph=true #To convert MACS14 XLS to Bedgraph and do Bedtools Coverage, set macs_14_bedgraph=true
#To process MNase-seq single-end data, use mnase_shift=true
#This will convert the bam file into a bed file using bedtools
#It will then extend each read to the fragment length estimated by Homer
#It will then reduce each read to the central 75bp portion and turn it back into a bam file
mnase_shift=false
#Step 10: Copy the files to the data directory
#To copy files to the data directory, set copy = true
#Once files are COPIED, set copy = false
copy_bam=true


##---------------------SET BIOPROJECT NUMBER----------------------------------------------
query=PRJNA146073

##---------------------SET GSE NUMBER HERE------------------------------------------------
#GSE number will identify the info and parallel command files being output in this script
des=GSE31221

##---------------------SET GENOME BUILD HERE----------------------------------------------
#Set genome build, mm10 or hg19
genome=mm10

#--------------------PREPROCESSING OPTIONS------------------------------------------------
#Set commands for adapter/quality trimming or quality masking here:
#For adapter/quality trimming: (set the proper adapter below)
atacadapter=false #Do adapter/quality trimming for ATAC seq datasets
chipadapter=true #Do adapter/quality trimming for ChIP seq datasets
rnaadapter=false #Do adapter/quality trimming for RNA seq datasets
mnaseadapter=false #Do adapter/quality trimming for MNase seq datasets

#Set adapter to be trimmed here
adapter="a GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAA"
#Set adapter to:
#illumina              for Illumina universal adapter
#nextera               for Nextera adapter
#small_rna             for Illumina Small RNA Adapter
#For a custom sequence, set adapter="a CUSTOMSEQUENCE"

#Set commands for quality masking instead of trimming (e.g., if low quality sequences are in the middle of the read)
atacquality=false #Do quality trimming for ATAC seq datasets
chipquality=false #Do quality trimming for ATAC seq datasets
rnaquality=false #Do quality trimming for ATAC seq datasets
mnasequality=false #Do quality trimming for ATAC seq datasets
#-----------------------------------------------------------------------------------------

#----------------------DUPLICATE REMOVAL OPTIONS------------------------------------------
#Set removed=true to remove duplicates (default)
#Set removed=false to mark duplicates, but keep the duplicates (e.g., for RNA-seq)
removed=true
#-----------------------------------------------------------------------------------------

#--------------Sets the working directories-----------------------------------------------
infodir=/mnt/data1/John/Pioneer-Factors/info
maindir=/mnt/data1/John/Pioneer_Factors
scriptdir=/mnt/data1/John/Pioneer-Factors
logdir=/mnt/data1/John/Pioneer-Factors/logs
filedir=/mnt/data1/John/Pioneer-Factors/sample_files
paralleldir=/mnt/data1/John/Pioneer-Factors/parallel_commands
r_dir=/mnt/data1/John/Pioneer-Factors/R_Scripts
datadir=/mnt/data1/VahediLab/PTF_Team/Data
datadir2=/mnt/data1/VahediLab/PTF_Team/Data2
dir0=$maindir
#-----------------------------------------------------------------------------------------

#--------------------CLEANS UP FILES------------------------------------------------------
rm -f $infodir/srr_files_${des}.txt
rm -f $paralleldir/srr_download_${des}.txt
rm -f $paralleldir/process_fastq_${des}.txt
rm -f $paralleldir/remove_duplicates_${des}.txt
rm -f $paralleldir/mark_duplicates_${des}.txt
rm -f $paralleldir/merge_bam_${des}.txt
rm -f $paralleldir/homer_key_file.txt
rm -f $paralleldir/macs_parallel_${des}.txt
#-----------------------------------------------------------------------------------------

#-------------Log File--------------------------------------------------------------------
start=`date +%s`
#rm -f $logdir/${des}_log.log
logfile=$logdir/${des}_log.log
exec > $logfile 2>&1
echo "-----------------------------------------------$start----------------------------------------------------------------" >> $logdir/${des}_log.log 
#------------------------------------------------------------------------------------------

##---------------------DOWNLOAD BIOPROJECT SAMPLE DATA------------------------------------
esearch -db sra -query $query | efetch --format docsum | xtract -pattern DocumentSummary -element Title > $infodir/sample_files_${des}.txt


##--------------Notes about SAMPLE FILE Format--------------------------------------------
#The format for the sample file needs to be in the following format:
#gsm,cell,species,seq,(mark),sequencer,type,replicate
#gsm should look like: GSM101236
#cell should look like: DP
#species should look like: Mus musculus
#seq should look like: ChIP_seq
#If seq is histone modification ChIP-seq, seq should look like ChIP_seq,H3K4me1
#If seq is transcription factor ChIP-seq, seq should look like ChIP_seq,Tcf1
#If seq is input DNA ChIP-seq, seq should look like ChIP_seq,Input
#investigator should be the name of the PI
#Note: Use unique identifier (such as investigator name) for input DNA to distinguish which experiment it belongs to
#sequencer should be either ABI or Illumina
#type should be either Paired or Single
#If there are biological replicates, replicate should look like Rep_1 or Rep_2
#If there are not biological replicates, replicate should look liks Rep_0

##---------------------CLEAN UP SAMPLE FILE-----------------------------------------------
#Use SED to clean up and properly format the SAMPLE FILE
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
sed -i -e 's/GSM773994: TCF7_ChIPSeq/GSM773994,EML,Mus musculus,ChIP_seq,Tcf1,Wu,Illumina,Single,Rep_0/g' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM773995: RUNX1_ChIPSeq_Rep1/GSM773995,EML,Mus musculus,ChIP_seq,Runx1,Wu,Illumina,Single,Rep_1/g' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM773996: RUNX1_ChIPSeq_Rep2/GSM773996,EML,Mus musculus,ChIP_seq,Runx1,Wu,Illumina,Single,Rep_2/g' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM773997: input DNA for TCF7/GSM773997,EML,Mus musculus,ChIP_seq,Input,Tcf1,Illumina,Single,Rep_0/g' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM773998: input DNA for RUNX1/GSM773998,EML,Mus musculus,ChIP_seq,Input,Runx1,Illumina,Single,Rep_0/g' $infodir/sample_files_${des}.txt


echo "#--------------------SAMPLE FILE-----------------------"
sed '' $infodir/sample_files_${des}.txt
echo "#--------------------------------------------------------"
echo "Notes about formatting:"
echo "The format for the sample file needs to be in the following format:"
echo "gsm,cell,species,seq,(mark),investigator,sequencer,type,replicate"
echo ""
echo "If the SAMPLE FILE is properly formated, set test1 = true"

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
#SRR123456.bam

#STAR generates genome index mm10
function star_generate_mm10 {
STAR --runMode genomeGenerate --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --genomeFastaFiles /mnt/data0/John/GRCm38p4_mm10_pa_only.fa --sjdbGTFfile /mnt/data0/John/gencode.vM6.annotation.gtf
}
#STAR generates genome index hg19
function star_generate_hg19 {
STAR --runMode genomeGenerate --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --genomeFastaFiles /mnt/data0/John/GRCh37.p13.genome.fa --sjdbGTFfile /mnt/data0/John/gencode.v19.annotation.gtf
}
#STAR alignment for ATAC-seq
function star_atac {
echo "STAR VERSION: `STAR --version`" >> $logdir/download_log_${des}.txt
if [[ "$genome" = mm* ]]; then
  echo "STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate" >> $logdir/download_log_${des}.txt
  STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate 
elif [[ "$genome" = hg* ]]; then
  echo "STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate" >> $logdir/download_log_${des}.txt
  STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate 
fi
mv Aligned.sortedByCoord.out.bam ${srr}.bam #Keeps the bam file that STAR generates
mv Log.final.out ${srr}_Log.final.out #Keeps a copy of the aligner log
}
#ABI alignment using BOWTIE
function abi-align {
echo "BOWTIE VERSION: `bowtie --version`" >> $logdir/download_log_${des}.txt
if [[ "$genome" = mm* ]]; then
  echo "bowtie -p 38 -S -C /mnt/data0/John/bowtie_mm10_colorspace/mm10_colorspace -f ${srr}_F3.csfasta -Q ${srr}_F3_QV.qual ${srr}.sam" >> $logdir/download_log_${des}.txt
  bowtie -p 38 -S -C /mnt/data0/John/bowtie_mm10_colorspace/mm10_colorspace -f ${srr}_F3.csfasta -Q ${srr}_F3_QV.qual ${srr}.sam
elif [[ "$genome" = hg* ]]; then
  echo "bowtie -p 38 -S -C /mnt/data0/John/bowtie_hg19_colorspace/hg19_c -f ${srr}_F3.csfasta -Q ${srr}_F3_QV.qual ${srr}.sam" >> $logdir/download_log_${des}.txt
  bowtie -p 38 -S -C /mnt/data0/John/bowtie_hg19_colorspace/hg19_c -f ${srr}_F3.csfasta -Q ${srr}_F3_QV.qual ${srr}.sam
fi

#cp /mnt/data0/John/bowtie_mm10_colorspace/mm10_colorspace.*.ebwt $wd #Copy ColorSpace Index file to working directory
mkdir -p tmp
echo "SAMTOOLS VERSION: `samtools --version`" >> $logdir/download_log_${des}.txt
samtools view -b -h ${srr}.sam > ${srr}_all_reads.bam #Convert sam to bam to save space
rm -f ${srr}.sam
echo "samtools view -b -h -F 4 ${srr}.sam > ${srr}_all_reads.bam #Keep only aligned reads" >> $logdir/download_log_${des}.txt
samtools view -b -h -F 4 ${srr}_all_reads.bam > ${srr}.bam #Keep only aligned reads
mv ${srr}.bam ${srr}_sort.bam
echo "------------SORTING BAM FILE-----------------" #Sort bam file according to coordinates
echo "PICARD VERSION: `java -jar $PICARD SortSam --version`" >> $logdir/download_log_${des}.txt
echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD SortSam INPUT=${srr}_sort.bam OUTPUT=${srr}.bam SORT_ORDER=coordinate" >> $logdir/download_log_${des}.txt
java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD SortSam INPUT=${srr}_sort.bam OUTPUT=${srr}.bam SORT_ORDER=coordinate
#cp ${srr}_F3.csfasta $datadir/$seq/reads
#cp ${srr}_F3_QV.qual $datadir/$seq/reads
}
#STAR alignment for ChIP-seq
function star_chip {
echo "STAR VERSION: `STAR --version`" >> $logdir/download_log_${des}.txt
if [[ "$genome" = mm* ]]; then
  echo "STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate" >> $logdir/download_log_${des}.txt
  STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate 
elif [[ "$genome" = hg* ]]; then
  echo "STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate" >> $logdir/download_log_${des}.txt
  STAR --runMode alignReads --alignIntronMax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 10000000000 --seedSearchStartLmax 30 --outFilterMultimapNmax 2 --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}.fastq --outSAMtype BAM SortedByCoordinate 
fi
mv Aligned.sortedByCoord.out.bam ${srr}.bam #Keeps the bam file that STAR generates
mv Log.final.out ${srr}_Log.final.out #Keeps a copy of the aligner log
}

#STAR alignment for RNA-seq 2-pass
function star_rna {
echo "STAR VERSION: `STAR --version`" >> $logdir/download_log_${des}.txt
if [[ "$genome" = mm* ]]; then
  echo "STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30" >> $logdir/download_log_${des}.txt
  STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/genome_GRCm38p4_M6 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30
elif [[ "$genome" = hg* ]]; then
  echo "STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30" >> $logdir/download_log_${des}.txt
  STAR --runMode alignReads --runThreadN 40 --genomeDir /mnt/data0/John/STAR_GRCh37p13_v19 --readFilesIn ${srr}_trimmed.fq --outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 30
fi
mv Aligned.sortedByCoord.out.bam ${srr}.bam #Keeps the bam file that STAR generates
mv Log.final.out ${srr}_Log.final.out #Keeps a copy of the aligner log
}
#STAR remove genome from memory
function star_remove {
STAR --genomeLoad Remove --genomeDir /mnt/data0/John/genome_GRCm38p4_M6
}

##---------------------SETTING FUNCTIONS FOR DUPLICATE REMOVAL----------------------------
#Duplicate removal should account for proper header information
#Generates commands for both marking and removal of duplicates

function mark_duplicates {
if [[ $(head -1 ${srr}.fastq | cut -c 1-4) = "@SRR" ]]; then
  #If file does not have header information regarding location of the read on the lane, will not look for optical duplicates
  echo "PICARD VERSION: `java -jar $PICARD MarkDuplicates --version`" >> $logdir/download_log_${des}.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_removed.bam READ_NAME_REGEX=null REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_${des}.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_marked.bam READ_NAME_REGEX=null METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_${des}.txt
else
  #If file has proper header information regarding location of the read on the lane, will look for optical duplicates
  echo "PICARD VERSION: `java -jar $PICARD MarkDuplicates --version`" >> $logdir/download_log_${des}.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_removed.bam REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_${des}.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_marked.bam METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_${des}.txt
fi
}

function mark_duplicates_abi {
  #Will not look for optical duplicates because ABI SOLiD 
  echo "PICARD VERSION: `java -jar $PICARD MarkDuplicates --version`" >> $logdir/download_log_${des}.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_removed.bam READ_NAME_REGEX=null REMOVE_DUPLICATES=true METRICS_FILE="`pwd`"/${srr}_metrics.txt" >> $paralleldir/remove_duplicates_${des}.txt
  echo "java -Xmx2g -jar $PICARD MarkDuplicates INPUT="`pwd`"/${srr}.bam OUTPUT="`pwd`"/${srr}_marked.bam READ_NAME_REGEX=null METRICS_FILE="`pwd`"/dupsMarked_${srr}_metrics.txt" >> $paralleldir/mark_duplicates_${des}.txt
}

##---------------------SETTING FUNCTIONS FOR PREPROCESSING--------------------------------
#Processes fastq files to remove adapters and low quality portions of reads
#Discards reads that are below 20bp
function adapter-trim {
  echo "TRIM_GALORE VERSION: `trim_galore --version`" >> $logdir/download_log_${des}.txt
  echo "trim_galore -o `pwd` --${adapter} `pwd`/${srr}.fastq" >> $logdir/download_log_${des}.txt
  cp ${srr}.fastq ${srr}_adapter.fastq
  trim_galore -o `pwd` --${adapter} `pwd`/${srr}.fastq
  mv ${srr}_trimmed.fq ${srr}.fastq
}

#Processes fastq files to mask low quality reads
function quality-trim {
  echo "Quality Trimming"
  echo "FASTX Toolkit (fastq_masker) version 0.0.14" >> $logdir/download_log_${des}.txt
  cp ${srr}.fastq ${srr}_quality.fastq
  echo "fastq_masker -q 20 -i `pwd`/${srr}.fastq -o `pwd`/${srr}_mask.fastq" >> $logdir/download_log_${des}.txt
  fastq_masker -q 20 -i `pwd`/${srr}.fastq -o `pwd`/${srr}_mask.fastq #Will mask low quality data
  mv ${srr}_mask.fastq ${srr}.fastq
}

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

#Combine fastq files for technical replicates
function fastq_combine {
allFastq=$(ls SRR*.fastq | grep -v STAR | grep -v adapter | grep -v marked | grep -v removed | grep -v 'sort' | grep -v reads | grep -v combine)
cat $allFastq > ${filename}.fastq
pigz -f -p 35 ${filename}.fastq
}
#Combines bam files and fastq files for biological replicates (Illumina)
function replicate-combine-illumina {
unset input
unset reads
for a in $(ls *GSM*.bam | grep -v STAR | grep -v adapter | grep -v marked | grep -v removed | grep -v 'sort' | grep -v reads | grep -v combine); do
  input=$(echo $input" INPUT="`pwd`/$a)
done
echo "PICARD VERSION MergeSamFiles: `java -jar $PICARD MergeSamFiles --version`" >> $logdir/download_log_${des}.txt
echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MergeSamFiles $input ""OUTPUT="`pwd`/${filename}_combine.bam >> $paralleldir/merge_bam_${des}.txt
allFastq=$(ls *GSM*.fastq.gz | grep -v STAR | grep -v adapter | grep -v marked | grep -v removed | grep -v 'sort' | grep -v reads)
cat $allFastq > ${filename}_combine.fastq.gz
}
#Combines bam files for biological replicates (ABI), but skips ABI fasta and qual files
function replicate-combine-abi {
unset input
unset reads
for a in $(ls *GSM*.bam | grep -v STAR | grep -v adapter | grep -v marked | grep -v removed | grep -v 'sort' | grep -v reads | grep -v combine); do
  input=$(echo $input" INPUT="`pwd`/$a)
done
echo "PICARD VERSION MergeSamFiles: `java -jar $PICARD MergeSamFiles --version`" >> $logdir/download_log_${des}.txt
echo "java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar $PICARD MergeSamFiles $input ""OUTPUT="`pwd`/${filename}_combine.bam >> $paralleldir/merge_bam_${des}.txt
}

#Processes single-end MNase data
function mnase-shift {
echo "bedtools bamtobed -i `pwd`/${filename}.bam > ${filename}_bam.bed" >> $logdir/download_log_${des}.txt 
bedtools bamtobed -i `pwd`/${filename}.bam > ${filename}_bam.bed
#---------------------Extending Reads-------------------------------
#Reads are extended by the estimated fragment length
echo "awk 'OFS="\t"{if ($6=="+") {start=$2; end=start+75} else if ($6=="-") {start=$3-75;end=$3} print $1,start,end,$4,$5,$6}' ${filename}_bam.bed > ${filename}_bam_extend.bed" >> $logdir/download_log_${des}.txt 
awk 'OFS="\t"{if ($6=="+") {start=$2; end=start+75} else if ($6=="-") {start=$3-75;end=$3} print $1,start,end,$4,$5,$6}' ${filename}_bam.bed > ${filename}_bam_extend.bed
#---------------------Shrinking Reads-------------------------------
#Reads are not reduced on each side until just the central 75bp portion is left
shiftLength=$(printf "%.0f" $(echo "scale=2;$fragLength/2" | bc))
echo "awk -v l=$shiftLength 'OFS="\t"{if ($6=="+") {start=$2+l; end=$3+l} else if ($6=="-") {start=$2-l;end=$3-l} print $1,start,end,$4,$5,$6}' ${filename}_bam_extend.bed > ${filename}_bam_extend_shift.bed" >> $logdir/download_log_${des}.txt 
awk -v l=$shiftLength 'OFS="\t"{if ($6=="+") {start=$2+l; end=$3+l} else if ($6=="-") {start=$2-l;end=$3-l} print $1,start,end,$4,$5,$6}' ${filename}_bam_extend.bed > ${filename}_bam_extend_shift.bed
#----------------------Convert to Bam-------------------------------
echo "sort -k1,1 -k2,2n ${filename}_bam_extend_shift.bed > ${filename}_bam_extend_shift_sort.bed" >> $logdir/download_log_${des}.txt 
sort -k1,1 -k2,2n ${filename}_bam_extend_shift.bed > ${filename}_bam_extend_shift_sort.bed
echo "bedToBam -i${filename}_bam_extend_shift_sort.bed -g $datadir/mm10.chrom.sizes > ${filename}_shift.bam" >> $logdir/download_log_${des}.txt 
bedToBam -i ${filename}_bam_extend_shift_sort.bed -g $datadir/mm10.chrom.sizes > ${filename}_shift.bam
}


#------------------------START OF SCRIPT--------------------------------------------------
if [[ $test1 = true ]]; then
echo ""; echo "-----SAMPLE FILE PROPERLY FORMATTED MOVING ONTO DOWNLOAD PREP-----"

##---------------------DOWNLOAD SRR NUMBERS ----------------------------------------------
#Will read in GSM numbers from SAMPLE FILE and download associated SRR numbers
#SRR File should have the following format:
#gsm,cell,species,seq,(mark),investigator,sequencer,type,replicate,wd:SRR1 SRR2 SRR3 ...
echo ""
echo "#--------------------SRR FILES---------------------------"
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
  cd $dir0/Data/ATAC_seq ; mkdir -p $cell/$investigator ; cd $cell/$investigator ; if [[ "$repnum" -ne 0 ]]; then mkdir -p $replicate; cd $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq at `pwd`"
  echo "$gsm,$cell,$species,$seq,$investigator,$sequencer,$type,$replicate,`pwd`:"`esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  investigator=$(echo $line | cut -d',' -f6)
  sequencer=$(echo $line | cut -d',' -f7)
  type=$(echo $line | cut -d',' -f8)
  replicate=$(echo $line | cut -d',' -f9)
  rep=$(echo $line | cut -d',' -f9 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f9| cut -d'_' -f2)
  cd $dir0/Data/ChIP_seq ; mkdir -p $cell/$mark/$investigator ; cd $cell/$mark/$investigator ; if [[ "$repnum" -ne 0 ]]; then mkdir -p $replicate; cd $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq $mark at `pwd`"
  echo "$gsm,$cell,$species,$seq,$mark,$investigator,$sequencer,$type,$replicate,`pwd`:"`esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; mkdir -p $cell/$investigator ; cd $cell/$investigator ; if [[ "$repnum" -ne 0 ]]; then mkdir -p $replicate; cd $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq at `pwd`"
  echo "$gsm,$cell,$species,$seq,$investigator,$sequencer,$type,$replicate,`pwd`:"`esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; mkdir -p $cell/$investigator ; cd $cell/$investigator ; if [[ "$repnum" -ne 0 ]]; then mkdir -p $replicate; cd $replicate; fi
  echo "Downloading SRR numbers for $gsm and making directory for $cell $seq at `pwd`"
  echo "$gsm,$cell,$species,$seq,$investigator,$sequencer,$type,$replicate,`pwd`:"`esearch -db sra -query "$gsm" | efetch --format runinfo | cut -d',' -f1 | grep SRR` >> $infodir/srr_files_${des}.txt
fi
done
echo "#--------------------------------------------------------"
echo "SRR File should have the following format gsm,cell,species,seq,(mark),investigator,sequencer,type,pwd:SRR1 SRR2 SRR3 ..." 
echo ""
sed '' $infodir/srr_files_${des}.txt
echo "If the SRR FILE is properly formated, set test2 = true"; echo ""

if [[ $test2 = true ]]; then
echo "----SETTING UP COMMANDS FOR DOWNLOAD-----"

###---------------------SRA TOOLKIT PARALLEL COMMANDS OUTPUT------------------------------
#Read in the SRR FILE, and set up parallel commands for later SRR download
echo "----------------FASTQ-DUMP---------------------"
while read line; do
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
wd=$(echo $line | cut -d',' -f9 | cut -d':' -f1)
if [[ "$seq" = ATAC* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "Download Paired-end Illumina data"
        fastq-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "Download Single-end Illumina data"
        fastq-download
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "Paired-end ABI data"
        abi-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "Download Single-end ABI data"
        abi-download
      fi
    fi
  done
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
  wd=$(echo $line | cut -d',' -f10 | cut -d':' -f1)
  cd $wd
  srr_line=$(echo $line | cut -d',' -f10 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "Download Paired-end Illumina data"
        fastq-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "Download Single-end Illumina data"
        fastq-download
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "Paired-end ABI data"
        abi-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "Download Single-end ABI data"
        abi-download
      fi
    fi
  done
fi
if [[ "$seq" = RNA* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "Download Paired-end Illumina data"
        fastq-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "Download Single-end Illumina data"
        fastq-download
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "Paired-end ABI data"
        abi-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "Download Single-end ABI data"
        abi-download
      fi
    fi
  done
fi
if [[ "$seq" = MNase* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = "Illumina" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "Download Paired-end Illumina data"
        fastq-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "Download Single-end Illumina data"
        fastq-download
      fi
    elif [[ "$sequencer" = "ABI" ]]; then
      if [[ "$type" = "Paired" ]]; then
        echo $srr "Paired-end ABI data"
        abi-download-paired
      elif [[ "$type" = "Single" ]]; then
        echo $srr "Download Single-end ABI data"
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
echo "Check to make sure download commands are correct before setting download = true"
echo ""

if [[ $download = true ]]; then
echo "--------CARRYING OUT THE ACTUAL DOWNLOAD----------"
echo "#--------------DOWNLOADING-------------------------" >> $logdir/download_log_${des}.txt
parallel --xapply --dryrun -j 30 -- < $paralleldir/srr_download_${des}.txt >> $logdir/download_log_${des}.txt
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
parallel --xapply --dryrun -j 30 -- < $paralleldir/srr_download_${des}.txt
parallel --xapply -j 30 -- < $paralleldir/srr_download_${des}.txt
echo ""; echo "-------------------------------------------------------"
echo "Recommended to stop here and check a few of the files for quality and adapter contamination using FASTQC"
echo "Set test3 = true when ready to continue"; echo ""
fi #download
                        
if [[ $test3 = true ]]; then
echo "MOVING ONTO PREPROCESSING"; echo ""
if [[ $process = true ]]; then
echo "--------------PREPROCESSING FILES--------------------"; echo ""
echo "#-------------PREPROCESSING------------------------" >> $logdir/download_log_${des}.txt
while read line; do
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
wd=$(echo $line | cut -d',' -f9 | cut -d':' -f1)
if [[ "$seq" = ATAC* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ $atacquality = true ]]; then
      echo $srr "Preprocessing: Quality Masking"    
      quality-trim
    fi
    if [[ $atacadapter = true ]]; then
      echo $srr "Preprocessing: Adapter and Quality Trimming"   
      adapter-trim
    fi
  done
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
  wd=$(echo $line | cut -d',' -f10 | cut -d':' -f1)
  cd $wd
  srr_line=$(echo $line | cut -d',' -f10 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ $chipadapter = true ]]; then
      echo $srr "Preprocessing: Quality Masking" 
      adapter-trim
    fi
    if [[ $chipquality = true ]]; then
      echo $srr "Preprocessing: Adapter and Quality Trimming" 
      quality-trim
    fi
  done
fi
if [[ "$seq" = RNA* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ $rnaquality = true ]]; then
      echo $srr "Preprocessing: Quality Masking" 
      quality-trim
    fi
    if [[ $rnaadapter = true ]]; then
      echo $srr "Preprocessing: Adapter and Quality Trimming" 
      adapter-trim
    fi
  done
fi
if [[ "$seq" = MNase* ]]
then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    echo $srr "quality or adapter trim"
    if [[ $mnasequality = true ]]; then
      echo $srr "Preprocessing: Quality Masking" 
      quality-trim
    fi
    if [[ $mnaseadapter = true ]]; then
      echo $srr "Preprocessing: Adapter and Quality Trimming" 
      adapter-trim
    fi
  done
fi
done < $infodir/srr_files_${des}.txt
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
echo "Once files are downloaded and preprocessed, set align = TRUE to align SRR files"
echo "-----------PREPROCESSING COMPLETED-------------"
fi #process

if [[ $align = true ]]; then
###---------------------ALIGNMENT---------------------------------
#Read in SRR files and align
echo ""; echo "--------------------ALIGNMENT--------------------"
echo "#--------------ALIGNMENT-------------------------" >> $logdir/download_log_${des}.txt
while read line; do
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
wd=$(echo $line | cut -d',' -f9 | cut -d':' -f1)
if [[ "$seq" = ATAC* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = Illumina ]]; then
      echo $srr "Align with STAR"
      star_atac
    elif [[ "$sequencer" = ABI ]]; then
      echo $srr "Align with BOWTIE"
      abi-align
    fi 
  done
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
  wd=$(echo $line | cut -d',' -f10 | cut -d':' -f1)
  cd $wd
  srr_line=$(echo $line | cut -d',' -f10 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = Illumina ]]; then
      echo $srr "Align with STAR"
      star_chip
    elif [[ "$sequencer" = ABI ]]; then
      echo $srr "Align with BOWTIE"
      abi-align
    fi 
  done
fi
if [[ "$seq" = RNA* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = Illumina ]]; then
      echo $srr "Align with STAR"
      star_rna
    elif [[ "$sequencer" = ABI ]]; then
      echo $srr "Align with BOWTIE"
      abi-align
    fi 
  done
fi
if [[ "$seq" = MNase* ]]; then
  cd $wd
  srr_line=$(echo $line | cut -d',' -f9 | cut -d':' -f2)
  IFS=' ' read -r -a srr_array <<< $srr_line
  for srr in "${srr_array[@]}"; do
    if [[ "$sequencer" = Illumina ]]; then
      echo $srr "Align with STAR"
      star_atac
    elif [[ "$sequencer" = ABI ]]; then
      echo $srr "Align with BOWTIE"
      abi-align
    fi 
  done
fi
done < $infodir/srr_files_${des}.txt
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
echo "#--------------------------------------------------------"
echo ""
echo ""; echo "Once files have been aligned, set duplicate = TRUE to mark and remove duplicates"; echo ""
fi #align

if [[ $duplicate = true ]]; then
echo "---------------REMOVE DUPLICATES--------------------"
while read line; do
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
wd=$(echo $line | cut -d',' -f9 | cut -d':' -f1)
if [[ "$seq" = ATAC* ]]; then
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
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  investigator=$(echo $line | cut -d',' -f6)
  sequencer=$(echo $line | cut -d',' -f7)
  type=$(echo $line | cut -d',' -f8)
  replicate=$(echo $line | cut -d',' -f9)
  rep=$(echo $line | cut -d',' -f9 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f9 | cut -d'_' -f2)
  wd=$(echo $line | cut -d',' -f10 | cut -d':' -f1)
  cd $wd
  srr_line=$(echo $line | cut -d',' -f10 | cut -d':' -f2)
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
if [[ "$seq" = MNase* ]]
then
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
done < $infodir/srr_files_${des}.txt

if [[ $removed = true ]]; then
  if [ -f $paralleldir/remove_duplicates_${des}.txt ]; then
  echo "-------REMOVING DUPLICATES------------"
  echo "#---------------REMOVING DUPLICATES-------------------------" >> $logdir/download_log_${des}.txt
  parallel --xapply --dryrun -j 30 -- < $paralleldir/remove_duplicates_${des}.txt >> $logdir/download_log_${des}.txt
  echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
  parallel --xapply --dryrun -j 30 -- < $paralleldir/remove_duplicates_${des}.txt
  parallel --xapply -j 30 -- < $paralleldir/remove_duplicates_${des}.txt
fi
else
if [ -f $paralleldir/mark_duplicates_${des}.txt ]; then
  echo "Marking, but not removing duplicates"
  echo "#-------------MARKING DUPLICATES-----------------------------" >> $logdir/download_log_${des}.txt
  parallel --xapply --dryrun -j 30 -- < $paralleldir/mark_duplicates_${des}.txt >> $logdir/download_log_${des}.txt
  echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
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
investigator=$(echo $line | cut -d',' -f5)
sequencer=$(echo $line | cut -d',' -f6)
type=$(echo $line | cut -d',' -f7)
replicate=$(echo $line | cut -d',' -f8)
rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
wd=$(echo $line | cut -d',' -f9 | cut -d':' -f1)
if [[ "$seq" = ATAC* ]]; then
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
if [[ "$seq" = ChIP* ]]; then
  seq=$(echo $line | cut -d',' -f4)
  mark=$(echo $line | cut -d',' -f5)
  investigator=$(echo $line | cut -d',' -f6)
  sequencer=$(echo $line | cut -d',' -f7)
  type=$(echo $line | cut -d',' -f8)
  replicate=$(echo $line | cut -d',' -f9)
  rep=$(echo $line | cut -d',' -f9 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f9 | cut -d'_' -f2)
  wd=$(echo $line | cut -d',' -f10 | cut -d':' -f1)
  cd $wd
  srr_line=$(echo $line | cut -d',' -f10 | cut -d':' -f2)
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
if [[ "$seq" = MNase* ]]
then
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
done < $infodir/srr_files_${des}.txt
echo "Now that duplicated sequences have been dealt with, set test4=true"
fi #duplicate

if [[ $test4 = true ]]; then
echo""; echo "COMBINING BAM FILES FOR TECHNICAL REPLICATES"
echo""; echo "---------------COMBINING BAM FILES--------------------"
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

if [[ $combine_bam = true ]]; then
echo "--------COMBINING BAM FILES--------"; echo ""
if [ -f $paralleldir/merge_bam_${des}.txt ]; then
echo "#----------------COMBINING TECHNICAL REPLICATES-----------------------" >> $logdir/download_log_${des}.txt
parallel --xapply --dryrun -j 30 -- < $paralleldir/merge_bam_${des}.txt >> $logdir/download_log_${des}.txt
parallel --xapply --dryrun -j 30 -- < $paralleldir/merge_bam_${des}.txt
parallel --xapply -j 30 -- < $paralleldir/merge_bam_${des}.txt
fi
fi #combine_bam

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
investigator=$(echo $line | cut -d',' -f5)
sequencer=$(echo $line | cut -d',' -f6)
type=$(echo $line | cut -d',' -f7)
replicate=$(echo $line | cut -d',' -f8)
rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell/$investigator ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  if [[ "$sequencer" = Illumina ]]; then
    echo "Combining fastq files for $filename" >> $logdir/download_log_${des}.txt
    fastq_combine
  elif [[ "$sequencer" = ABI ]]; then
    echo "Skipping ${filename} fastq combine because ABI" >> $logdir/download_log_${des}.txt
  fi 
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
  if [[ "$sequencer" = Illumina ]]; then
    echo "Combining fastq files for $filename" >> $logdir/download_log_${des}.txt
    fastq_combine
  elif [[ "$sequencer" = ABI ]]; then
    echo "Skipping ${filename} fastq combine because ABI" >> $logdir/download_log_${des}.txt
  fi 
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  if [[ "$sequencer" = Illumina ]]; then
    echo "Combining fastq files for $filename" >> $logdir/download_log_${des}.txt
    fastq_combine
  elif [[ "$sequencer" = ABI ]]; then
    echo "Skipping ${filename} fastq combine because ABI" >> $logdir/download_log_${des}.txt
  fi 
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator ; if [[ "$repnum" -ne 0 ]]; then cd $replicate; fi
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  if [[ "$sequencer" = Illumina ]]; then
    echo "Combining fastq files for $filename" >> $logdir/download_log_${des}.txt
    fastq_combine
  elif [[ "$sequencer" = ABI ]]; then
    echo "Skipping ${filename} fastq combine because ABI" >> $logdir/download_log_${des}.txt
  fi 
fi
done
echo "#--------------------------------------------------------"
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
fi #combine_fastq

if [[ $replicate_copy = true ]]; then
echo "Combining Biological Replicates"; echo ""
echo "#--------------------Combining Biological Replicates---------------------------"
echo "#Part 1 Copying Files"
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
  cd $dir0/Data/ATAC_seq ; cd $cell/$investigator
  if [[ "$repnum" -ne 0 ]]; then 
    cd $replicate
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "Copying Rep $repnum fastq to upper directory"
    cp ${filename}.bam ..
    if [[ "$sequencer" = Illumina ]]; then
      cp ${filename}.fastq.gz ..
    elif [[ "$sequencer" = ABI ]]; then
      echo "Skipping ${filename} because ABI"
    fi 
  fi
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
  cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator
  if [[ "$repnum" -ne 0 ]]; then
    cd $replicate
    filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
    echo "Copying $cell $mark Rep $repnum fastq to upper directory"
    cp ${filename}.bam ..
    if [[ "$sequencer" = Illumina ]]; then
      cp ${filename}.fastq.gz ..
    elif [[ "$sequencer" = ABI ]]; then
      echo "Skipping ${filename} because ABI"
    fi 
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  if [[ "$repnum" -ne 0 ]]; then
    cd $replicate
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "Copying Rep $repnum fastq to upper directory"
    cp ${filename}.bam ..
    if [[ "$sequencer" = Illumina ]]; then
      cp ${filename}.fastq.gz ..
    elif [[ "$sequencer" = ABI ]]; then
      echo "Skipping ${filename} because ABI"
    fi 
  fi
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  if [[ "$repnum" -ne 0 ]]; then
    cd $replicate
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "Copying Rep $repnum fastq to upper directory"
    cp ${filename}.bam ..
    if [[ "$sequencer" = Illumina ]]; then
      cp ${filename}.fastq.gz ..
    elif [[ "$sequencer" = ABI ]]; then
      echo "Skipping ${filename} because ABI"
    fi 
  fi
fi
done
fi #replicate_copy

rm -f $paralleldir/merge_bam_${des}.txt
rm -f $infodir/sample_files_${des}_noreps.txt

if [[ $test5 = true ]]; then
echo "#----------Making New Sample File and Generating Homer Key File---------------------------"
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
  cd $dir0/Data/ATAC_seq ; cd $cell/$investigator
  if [[ "$repnum" = 0 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
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
  cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator
  if [[ "$repnum" = 0 ]]; then 
    filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    if echo "$mark" | grep -q Input; then
      echo "Input found, not making Homer Tag Directory: `pwd`/${filename}.bam"
    else
      echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
    fi
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    if echo "$mark" | grep -q Input; then
      echo "Input found, not making Homer Tag Directory: `pwd`/${filename}.bam"
    else
      echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
    fi
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  if [[ "$repnum" = 0 ]]; then 
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  echo "$line" >> $infodir/sample_files_${des}_noreps.txt
  echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  if [[ "$repnum" = 0 ]]; then 
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  echo "$line" >> $infodir/sample_files_${des}_noreps.txt
  echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    echo "$line" >> $infodir/sample_files_${des}_noreps.txt
    echo -e "${filename}"'\t'`pwd`"/${filename}.bam" >> $paralleldir/homer_key_file.txt
  fi
fi
done
sed '' $infodir/sample_files_${des}_noreps.txt
echo "#------------------------------------------------------------------"
echo "To combine biological replicates, set replicate_combine=true"
echo ""

if [[ $replicate_combine = true ]]; then
echo "#-------------COMBINING BIOLOGICAL REPLICATES---------------------------"
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
  cd $dir0/Data/ATAC_seq ; cd $cell/$investigator
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    if [[ "$sequencer" = Illumina ]]; then
      replicate-combine-illumina
    elif [[ "$sequencer" = ABI ]]; then
      replicate-combine-abi
    fi 
  fi
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
  cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
    if [[ "$sequencer" = Illumina ]]; then
      replicate-combine-illumina
    elif [[ "$sequencer" = ABI ]]; then
      replicate-combine-abi
    fi 
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    if [[ "$sequencer" = Illumina ]]; then
      replicate-combine-illumina
    elif [[ "$sequencer" = ABI ]]; then
      replicate-combine-abi
    fi 
  fi
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  if [[ "$repnum" = 1 ]]; then 
    filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
    if [[ "$sequencer" = Illumina ]]; then
      replicate-combine-illumina
    elif [[ "$sequencer" = ABI ]]; then
      replicate-combine-abi
    fi 
  fi
fi
done
echo "#--------------------------------------------------------"
echo ""
sed '' $infodir/sample_files_${des}_noreps.txt

echo "COMBINING BAM FILES"; echo ""
if [ -f $paralleldir/merge_bam_${des}.txt ]; then
echo "#---------------COMBINING BIO REPS--------------------------" >> $logdir/download_log_${des}.txt
parallel --xapply --dryrun -j 30 -- < $paralleldir/merge_bam_${des}.txt >> $logdir/download_log_${des}.txt
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
parallel --xapply --dryrun -j 30 -- < $paralleldir/merge_bam_${des}.txt
parallel --xapply -j 30 -- < $paralleldir/merge_bam_${des}.txt
echo "#-------------Clean Up Directories---------------------------"
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
  cd $dir0/Data/ATAC_seq ; cd $cell
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  if [[ "$repnum" -ne 0 ]]; then
    if [[ "$repnum" -eq 1 ]]; then 
      if [[ "$sequencer" = Illumina ]]; then
        rm -f ${filename}.bam
        rm -f ${filename}.fastq.gz
        mv ${filename}_combine.bam ${filename}.bam
        mv ${filename}_combine.fastq.gz ${filename}.fastq.gz
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
        mv ${filename}_combine.bam ${filename}.bam
      fi 
    else
      if [[ "$sequencer" = Illumina ]]; then
        rm -f ${filename}.bam
        rm -f ${filename}.fastq.gz
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
      fi 
    fi  
  fi
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
  cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator
  filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
  if [[ "$repnum" -ne 0 ]]; then
    if [[ "$repnum" -eq 1 ]]; then 
      if [[ "$sequencer" = Illumina ]]; then
        rm -f ${filename}.bam
        rm -f ${filename}.fastq.gz
        mv ${filename}_combine.bam ${filename}.bam
        mv ${filename}_combine.fastq.gz ${filename}.fastq.gz
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
        mv ${filename}_combine.bam ${filename}.bam
      fi 
    else
      if [[ "$sequencer" = Illumina ]]; then
        rm -f ${filename}.bam
        rm -f ${filename}.fastq.gz
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
      fi 
    fi  
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  if [[ "$repnum" -ne 0 ]]; then
    if [[ "$repnum" -eq 1 ]]; then 
      if [[ "$sequencer" = Illumina ]]; then
        rm -f ${filename}.bam
        rm -f ${filename}.fastq.gz
        mv ${filename}_combine.bam ${filename}.bam
        mv ${filename}_combine.fastq.gz ${filename}.fastq.gz
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
        mv ${filename}_combine.bam ${filename}.bam
      fi 
    else
      if [[ "$sequencer" = Illumina ]]; then
        rm -f ${filename}.bam
        rm -f ${filename}.fastq.gz
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
      fi 
    fi  
  fi
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  if [[ "$repnum" -ne 0 ]]; then
    if [[ "$repnum" -eq 1 ]]; then 
      if [[ "$sequencer" = Illumina ]]; then
        rm -f ${filename}.bam
        rm -f ${filename}.fastq.gz
        mv ${filename}_combine.bam ${filename}.bam
        mv ${filename}_combine.fastq.gz ${filename}.fastq.gz
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
        mv ${filename}_combine.bam ${filename}.bam
      fi 
    else
      if [[ "$sequencer" = Illumina ]]; then
        rm -f ${filename}.bam
        rm -f ${filename}.fastq.gz
      elif [[ "$sequencer" = ABI ]]; then
        rm -f ${filename}.bam
      fi 
    fi  
  fi
fi
done
fi
fi #replicate_combine

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
lines=$(wc -l $infodir/sample_files_${des}_noreps.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}_noreps.txt)
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
  cd $dir0/Data/ATAC_seq ; cd $cell ; mkdir -p macs
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/$des/$filename/tagInfo.txt" | cut -d"=" -f2 | xargs)
  echo "`pwd`/${filename}.bam `pwd`/macs/${filename}" "$fragLength">> $paralleldir/macs_parallel_${des}.txt
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
  filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
  if echo "$mark" | grep -q Input; then
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator
    echo "Input found: `pwd`/${filename}.bam"
  else
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator ; mkdir -p macs
    fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/$des/$filename/tagInfo.txt" | cut -d"=" -f2 | xargs)
    echo "`pwd`/${filename}.bam `pwd`/macs/${filename}" "$fragLength">> $paralleldir/macs_parallel_${des}.txt
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
done

rm -f $paralleldir/macs_parallel_chip_${des}.txt
rm -f $paralleldir/macs_parallel_atac_${des}.txt

##---------------------CLEAN UP MACS14 PARALLEL COMMANDS-----------------------------------------------
#Here are sample sed commands:
#Replace text with NewText
###sed -i -e 's/text/NewText/g' $paralleldir/macs_parallel_${des}.txt
#Delete lines with pattern
###sed -i -e '/pattern/d' $paralleldir/macs_parallel_${des}.txt
#Print lines matching pattern
###sed -i -n '/pattern;/p' $paralleldir/macs_parallel_${des}.txt
sed ' ' $paralleldir/macs_parallel_${des}.txt > $paralleldir/macs_parallel_chip_${des}.txt
sed -i -e 's, 118, 118 /mnt/data1/John/Pioneer_Factors/Data/ChIP_seq/EML/Input/Runx1/EML_ChIP_seq_Input_Runx1_mm10_GSM773998.bam,g' $paralleldir/macs_parallel_chip_${des}.txt
sed -i -e 's, 98, 98 /mnt/data1/John/Pioneer_Factors/Data/ChIP_seq/EML/Input/Tcf1/EML_ChIP_seq_Input_Tcf1_mm10_GSM773997.bam,g' $paralleldir/macs_parallel_chip_${des}.txt

#sed -i -n '/pattern;/p' $paralleldir/macs_parallel_${des}.txt

echo ""
echo "#------------------CHECK MACS14 PARALLEL COMMANDSf------------------------"
sed '' $paralleldir/macs_parallel_${des}.txt
echo "#--------------------------------------------------------"
echo ""
echo "#--------------------------------------------------------"
echo ""
echo "-----------------CHECK MACS14 ATAC Seq PARALLEL COMMANDS-----------------"
sed '' $paralleldir/macs_parallel_atac_${des}.txt
echo "#--------------------------------------------------------"
echo "-----------------CHECK MACS14 ChIP Seq PARALLEL COMMANDS-----------------"
sed '' $paralleldir/macs_parallel_chip_${des}.txt
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
if [ -f $paralleldir/macs_parallel_chip_${des}.txt ]; then
echo "#-----------------MACS PEAK CALLING--------------------------" >> $logdir/download_log_${des}.txt
parallel --dryrun --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip_${des}.txt "macs14 -t {1} -n {2} -c {4} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"  >> $logdir/download_log_${des}.txt
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
parallel --dryrun --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip_${des}.txt "macs14 -t {1} -n {2} -c {4} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"
parallel --xapply -j 38 --colsep ' ' -a $paralleldir/macs_parallel_chip_${des}.txt "macs14 -t {1} -n {2} -c {4} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"
fi
fi #macs_14_chip 

##----------------------MACS Peak Calling ATAC seq----------------------------------------
if [[ $macs_14_atac = true ]]; then
echo ""; echo "--------------MACS Peak Calling ATAC seq---------------"
if [ -f $paralleldir/macs_parallel_atac_${des}.txt ]; then
echo "#-------------------MACS PEAK CALLING--------------------------" >> $logdir/download_log_${des}.txt
parallel --xapply --dryrun -j 35 --colsep ' ' -a $paralleldir/macs_parallel_atac_${des}.txt "macs14 -t {1} -n {2} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1" >> $logdir/download_log_${des}.txt
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
parallel --xapply --dryrun -j 35 --colsep ' ' -a $paralleldir/macs_parallel_atac_${des}.txt "macs14 -t {1} -n {2} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"
parallel --xapply -j 35 --colsep ' ' -a $paralleldir/macs_parallel_atac_${des}.txt "macs14 -t {1} -n {2} --bw {3} -f BAM -g mm -p $MACSpvalue -w --single-profile >> {2}.log 2>&1"
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
investigator=$(echo $line | cut -d',' -f5)
sequencer=$(echo $line | cut -d',' -f6)
type=$(echo $line | cut -d',' -f7)
replicate=$(echo $line | cut -d',' -f8)
rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell ; cd macs
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
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
  filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
  if echo "$mark" | grep -q Input; then
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator
    echo "Skipping Input"
  else
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator ; cd macs
    for a in `ls *model.r`; do
      echo "Generating MACS model figure PDF for ${filename}"
      Rscript --verbose $a
    done
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  echo "No Peak Calling for RNA-seq"
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  echo "No Peak Calling for MNase-seq"
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
done
fi #macs_14_model

#---------------------CHECK MACS ERRORS---------------------------------------------------
#Not in this current version

#---------------------Normalize WIG FILE--------------------------------------------------
if [[ $macs_14_bigwig = true ]]; then
echo ""; echo "-------------Normalize WIG FILE---------------"
echo "#----------------NORMALIZING WIG FILE--------------------------" >> $logdir/download_log_${des}.txt
lines=$(wc -l $infodir/sample_files_${des}_noreps.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}_noreps.txt)
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
  cd $dir0/Data/ATAC_seq ; cd $cell ; cd macs
  filename=${cell}_${seq}_${investigator}_${genome}_${gsm}
  wd=`pwd`
  if [ ! -f ${filename}_normalized.wig ]; then
  if [ -f ${filename}_treat_afterfiting_all.wig.gz ]; then
    gzip -d ${filename}_treat_afterfiting_all.wig.gz
    readsTotal=$(samtools view -c ../${filename}.bam)
    reads=$(bc <<< "scale=2;$readsTotal/1000000")
    echo "${filename} normalizing of wig file, reads in millions: $reads" >> $logdir/download_log_${des}.txt
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
  investigator=$(echo $line | cut -d',' -f6)
  sequencer=$(echo $line | cut -d',' -f7)
  type=$(echo $line | cut -d',' -f8)
  replicate=$(echo $line | cut -d',' -f9)
  rep=$(echo $line | cut -d',' -f9 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f9 | cut -d'_' -f2)
  filename=${cell}_${seq}_${mark}_${investigator}_${genome}_${gsm}
  if echo "$mark" | grep -q Input; then
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator
    echo "Skipping Input"
  else
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator ; cd macs
      wd=`pwd`
      if [ ! -f ${filename}_normalized.wig ]; then
      if [ -f ${filename}_treat_afterfiting_all.wig.gz ]; then
        gzip -d ${filename}_treat_afterfiting_all.wig.gz
        readsTotal=$(samtools view -c ../${filename}.bam)
        reads=$(bc <<< "scale=2;$readsTotal/1000000")
        echo "${filename} normalizing of wig file, reads in millions: $reads" >> $logdir/download_log_${des}.txt
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
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  echo "No Peak Calling for RNA-seq"
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  echo "No Peak Calling for MNase-seq"
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
done
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
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
investigator=$(echo $line | cut -d',' -f5)
sequencer=$(echo $line | cut -d',' -f6)
type=$(echo $line | cut -d',' -f7)
replicate=$(echo $line | cut -d',' -f8)
rep=$(echo $line | cut -d',' -f8 | cut -d'_' -f1)
repnum=$(echo $line | cut -d',' -f8 | cut -d'_' -f2)
if [[ "$seq" = ATAC* ]]; then
  cd $dir0/Data/ATAC_seq ; cd $cell ; cd macs
  filename=${cell}_${seq}_${investigator}_mm10_${gsm}
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
  investigator=$(echo $line | cut -d',' -f6)
  sequencer=$(echo $line | cut -d',' -f7)
  type=$(echo $line | cut -d',' -f8)
  replicate=$(echo $line | cut -d',' -f9)
  rep=$(echo $line | cut -d',' -f9 | cut -d'_' -f1)
  repnum=$(echo $line | cut -d',' -f9 | cut -d'_' -f2)
  filename=${cell}_${seq}_${mark}_${investigator}_mm10_${gsm}
  if echo "$mark" | grep -q Input; then
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator
    echo "Skipping Input"
  else
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator ; cd macs
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
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  echo "No Peak Calling for RNA-seq"
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  echo "No Peak Calling for MNase-seq"
  #---------NO PEAK CALLING FOR RNA-SEQ------------------------
fi
done
fi #macs_14_bedgraph

fi #macs_14

#----------------------Processing MNase_seq bam files-----------------------------------------
if [[ $mnase_shift = true ]]; then
echo "#----------------Generating Shifted Bam file for MNase------------------" >> $logdir/download_log_${des}.txt
echo ""; echo "#--------------------Generating Shifted Bam file for MNase---------------------------"
lines=$(wc -l $infodir/sample_files_${des}_noreps.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}_noreps.txt)
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
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  filename=${cell}_${seq}_${investigator}_mm10_${gsm}
  fragLength=$(sed -n '3p' < "$dir0/Analysis/Homer/Tag_Directories/$des/$filename/tagInfo.txt" | cut -d"=" -f2 | xargs)
  mnase-shift
fi
done
echo "#--------------------------------------------------------"
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
echo ""; echo "#--------------------Copying Shifted Bam to Data Dir---------------------------"
echo "#-----------------Shifted Bam Copy------------------------" >> $logdir/download_log_${des}.txt 
lines=$(wc -l $infodir/sample_files_${des}_noreps.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}_noreps.txt)
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
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  filename=${cell}_${seq}_${investigator}_mm10_${gsm}
  DATE=`date +%Y-%m-%d:%H:%M:%S`
  echo "Copying Shifted Bam: ${filename}_shift.bam to Data Dir: $datadir2 on $DATE" >> $logdir/download_log_${des}.txt 
  samtools index ${filename}_shift.bam
  cp ${filename}_shift.bam $datadir2/MNase_seq/bam
  cp ${filename}_shift.bam.bai $datadir2/MNase_seq/bam
fi
done
echo "#--------------------------------------------------------" >> $logdir/download_log_${des}.txt
echo "#--------------------------------------------------------"
fi #mnase_shift

if [[ $copy_bam = true ]]; then
echo "#------------------COPYING TO DATA DIR------------------" >> $logdir/download_log_${des}.txt
echo "COPYING BAM FILES TO DATA DIR"; echo ""
echo "#--------------------Copy Bam File to Data Dir---------------------------"
lines=$(wc -l $infodir/sample_files_${des}_noreps.txt | cut -d' ' -f1)
for ((i=1; i<=$lines; i++)); do
line=$(sed -n "${i}p" < $infodir/sample_files_${des}_noreps.txt)
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
  cd $dir0/Data/ATAC_seq ; cd $cell
  filename=${cell}_${seq}_${investigator}_mm10_${gsm}
  echo "Copying Bam and MACS files for $filename to Data Directory $datadir2" >> $logdir/download_log_${des}.txt 
  samtools index ${filename}.bam
  cp ${filename}.bam $datadir2/ATAC_seq/bam
  cp ${filename}.bam.bai $datadir2/ATAC_seq/bam
  cp -R macs $datadir2/ATAC_seq/
  cp macs/${filename}_normalized.bw $datadir2/ATAC_seq/bw
  if [[ "$sequencer" = Illumina ]]; then
    cp ${filename}.fastq.gz $datadir2/ATAC_seq/fastq
  fi 
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
  cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator
  filename=${cell}_${seq}_${mark}_${investigator}_mm10_${gsm}
  if echo "$mark" | grep -q Input; then
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator
    echo "Skipping Input"
  else
    cd $dir0/Data/ChIP_seq ; cd $cell/$mark/$investigator
    echo "Copying Bam and MACS files for $filename to Data Directory $datadir2" >> $logdir/download_log_${des}.txt 
    samtools index ${filename}.bam
    cp ${filename}.bam $datadir2/ChIP_seq/bam
    cp ${filename}.bam.bai $datadir2/ChIP_seq/bam
    cp -R macs $datadir2/ChIP_seq/
    cp macs/${filename}_normalized.bw $datadir2/ChIP_seq/bw
    if [[ "$sequencer" = Illumina ]]; then
      cp ${filename}.fastq.gz $datadir2/ChIP_seq/fastq
    fi 
  fi
fi
if [[ "$seq" = RNA* ]]; then
  cd $dir0/Data/RNA_seq ; cd $cell/$investigator
  filename=${cell}_${seq}_${investigator}_mm10_${gsm}
  echo "Copying Bam and MACS files for $filename to Data Directory $datadir2" >> $logdir/download_log_${des}.txt 
  samtools index ${filename}.bam
  cp ${filename}.bam $datadir2/RNA_seq/bam
  cp ${filename}.bam.bai $datadir2/RNA_seq/bam
  if [[ "$sequencer" = Illumina ]]; then
    cp ${filename}.fastq.gz $datadir2/RNA_seq/fastq
  fi 
fi
if [[ "$seq" = MNase* ]]; then
  cd $dir0/Data/MNase_seq ; cd $cell/$investigator
  filename=${cell}_${seq}_${investigator}_mm10_${gsm}
  echo "Copying Bam and MACS files for $filename to Data Directory $datadir2" >> $logdir/download_log_${des}.txt 
  samtools index ${filename}.bam
  cp ${filename}.bam $datadir2/MNase_seq/bam
  cp ${filename}.bam.bai $datadir2/MNase_seq/bam
  if [[ "$sequencer" = Illumina ]]; then
    cp ${filename}.fastq.gz $datadir2/MNase_seq/fastq
  fi
fi
done
echo "#--------------------------------------------------------"
echo ""; echo "Finished copying"; echo ""
DATE=`date +%Y-%m-%d:%H:%M:%S`
echo "----------------------------------------------------$DATE---------------------------------------------------------------------------" >> $logdir/download_log_${des}.txt 
echo "Notes: Runx1 and Tcf1 ChIP-seq in EML (progenitor) cell line. Tcf1 dataset is good quality, Runx1 dataset is lower quality, although trimmed." >> $logdir/download_log_${des}.txt 
echo "-------------------------------------------------------------------------------------------------------------------------------------" >> $logdir/download_log_${des}.txt  
cp $logdir/download_log_${des}.txt $datadir2/logs
fi #copy_bam

#perl -ne 'print unless $seen{$_}++' $logdir/download_log_${des}.txt

fi #test5
fi #test4
fi #test3
fi #test2
fi #test1

