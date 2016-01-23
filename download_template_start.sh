#!/bin/bash
#This script will download the sample and SRR numbers for a given bioproject, and do the following:
#1. Process and align the files
#2. Combine technical and biological replicates
#3. Do peak calling or other processing
#4. Generate normalized bigwig files
#5. Move files to the data directory

#-------------------------HOW TO and OPTIONS----------------------------------------------
#Step 1: Add the BIOPROJECT NUMBER of the sample to be downloaded
#Step 2: Add the GSE NUMBER of the sample to be downloaded 
#Step 3: Run the script, use SED to properly format the file according to guidelines
#Step 4: Check to make sure the SRR FILE is properly formatted
#Step 5: Use FASTQC to check the fastq files for quality and adapter contaminations
#Step 6: Set the options for preprocessing (below) and then set test3=true to continue
#Step 7: When files have been processed, aligned, deduplicated, and technical replicates are now ready to be combined, set test4=true
#Step 8: To combine biological replicates, set replicate_copy=true to copy the files for the replicates into the right place
#Step 9: Choose or omit Peak Calling options below
#Step 10: Copy the files to the data directory

test1=true #Once the SAMPLE FILE is properly FORMATTED, set test1=true
test2=true #Once the SRR FILE is properly FORMATTED, set test2=true
  download=false #To CARRY OUT the download, set download=true. Once DOWNLOADED, set download=false to avoid redownloading it.
test3=true #Once you have checked the files with FASTQC, set test3=true
  process=false #To CARRY OUT the preprocessing, set process=true. Once files are PROCESSED, set process=false
  align=false #To ALIGN, set align = true. Once ALIGNED, set align=false
  duplicate=false #To remove DUPLICATES, set duplicate=true. Once DUPLICATES have been removed, set duplicate=false. #Note: If you want to keep duplicates, set removed=false (below) before running duplicate=true
test4=true
  combine_bam=true #To CARRY OUT the actual combining of bam files, set combine_bam=true. Once files are COMBINED, set combine_bam=false
  combine_fastq=true #To COMBINE fastq files, set combine_fastq=true. If files already COMBINED, set combine_fastq=false
test5=true #Once bam and fastq files combined for technical replicates are combined, set test5=true
bio_replicate=true #To COMBINE biological replicates, set bio_replicate=true. Set to false if want to keep biological replicates separate.
  replicate_copy=true #To copy biological replicates to the right place, set replicate_copy=true
  replicate_combine=true #To combine bam biological replicate files, set replicate_combine=false
homer=false #To CARRYOUT Homer makeTagDirectory, set homer=true. Once tag directories are CREATED, set homer=false
macs_14=false #For MACS14 Peak Calling, set macs_14=true. #Check to make sure the MACS14 commands file is correct (below) before setting an option and CARRYING OUT the peak calling
  macs_14_chip=false #For TF ChIP-seq (with Input) MACS14 set macs_14_chip=true
  macs_14_atac=false #For ATAC-seq MACS14 for single-end data, set macs_14_chip=true
  macs_14_model=false #To build Model PDF from R script, set macs_14_model=true
  macs_14_bigwig=false #To normalize BigWig file from MACS, set macs_14_bigwig=true
  macs_14_bedgraph=false #To convert MACS14 XLS to Bedgraph and do Bedtools Coverage, set macs_14_bedgraph=true
macs_2=false #For MACS2 Peak Calling, set macs_2=true. Check to make sure the MACS2 commands file is correct (below) before setting an option and CARRYING OUT the peak calling
  macs_2_chip=false #For TF ChIP-seq (with Input) MACS2 set macs_2_chip=true
  macs_2_atac=false #For ATAC-seq MACS2 for paired-end data, set macs_2_chip=true
  macs_2_model=false #To build Model PDF from R script, set macs_2_model=true
  macs_2_bigwig=false #To normalize BigWig file from MACS, set macs_14_bigwig=true
  macs_2_bedgraph=false #To convert MACS14 XLS to Bedgraph and do Bedtools Coverage, set macs_14_bedgraph=true
#This will convert the bam file into a bed file using bedtools
#It will then extend each read to the fragment length estimated by Homer
#It will then reduce each read to the central 75bp portion and turn it back into a bam file
mnase_shift=false #To process MNase-seq single-end data, use mnase_shift=true
#To copy files to the data directory, set copy = true
#Once files are COPIED, set copy = false
copy_bam=false

echo "#--------------SET BIOPROJECT NUMBER-----------------------------------"
query=PRJNA306754
echo "Bioproject number set as: $query"

echo "#--------------SET GSE NUMBER -----------------------------------------"
#GSE number will identify the info and parallel command files being output in this script
des=GSE76268
echo "GSE number set as: $des"

echo "#--------------SET GENOME BUILD----------------------------------------"
#Set genome build: mm10 or hg19
genome=hg19
echo "Genome number set as: $genome"

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
start=`date +%Y-%m-%d:%H:%M:%S`
#rm -f $logdir/${des}_log.log
logfile=$logdir/${des}_log.log
exec > $logfile 2>&1
echo "---------------------------------------------Start of Script: $start----------------------------------------------------------------" >> $logdir/${des}_log.log 
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
sed -i -e '/RNA-seq/d' $infodir/sample_files_${des}.txt
sed -i -e '/Acinar/d' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM1978248: ATAC-seq Beta3; Homo sapiens; OTHER/GSM1978248,Beta,Homo sapiens,ATAC_seq,Wang,Illumina,Paired,Rep_1/g' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM1978247: ATAC-seq Beta2; Homo sapiens; OTHER/GSM1978247,Beta,Homo sapiens,ATAC_seq,Wang,Illumina,Paired,Rep_2/g' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM1978246: ATAC-seq Beta1; Homo sapiens; OTHER/GSM1978246,Beta,Homo sapiens,ATAC_seq,Wang,Illumina,Paired,Rep_3/g' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM1978245: ATAC-seq Alpha3; Homo sapiens; OTHER/GSM1978245,Alpha,Homo sapiens,ATAC_seq,Wang,Illumina,Paired,Rep_1/g' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM1978244: ATAC-seq Alpha2; Homo sapiens; OTHER/GSM1978244,Alpha,Homo sapiens,ATAC_seq,Wang,Illumina,Paired,Rep_2/g' $infodir/sample_files_${des}.txt
sed -i -e 's/GSM1978243: ATAC-seq Alpha1; Homo sapiens; OTHER/GSM1978243,Alpha,Homo sapiens,ATAC_seq,Wang,Illumina,Paired,Rep_3/g' $infodir/sample_files_${des}.txt
#sed -i -e 's/GSM773995: RUNX1_ChIPSeq_Rep1/GSM773995,EML,Mus musculus,ChIP_seq,Runx1,Wu,Illumina,Single,Rep_1/g' $infodir/sample_files_${des}.txt
#sed -i -e 's/GSM773996: RUNX1_ChIPSeq_Rep2/GSM773996,EML,Mus musculus,ChIP_seq,Runx1,Wu,Illumina,Single,Rep_2/g' $infodir/sample_files_${des}.txt
#sed -i -e 's/GSM773997: input DNA for TCF7/GSM773997,EML,Mus musculus,ChIP_seq,Input,Tcf1,Illumina,Single,Rep_0/g' $infodir/sample_files_${des}.txt
#sed -i -e 's/GSM773998: input DNA for RUNX1/GSM773998,EML,Mus musculus,ChIP_seq,Input,Runx1,Illumina,Single,Rep_0/g' $infodir/sample_files_${des}.txt


echo "#--------------SAMPLE FILE-----------------------"
sed '' $infodir/sample_files_${des}.txt
echo "#--------------------------------------------------------"
echo "Notes about formatting:"
echo "The format for the sample file needs to be in the following format:"
echo "gsm,cell,species,seq,(mark),investigator,sequencer,type,replicate"

if [[ $test1 = true ]]; then
  . ${scriptdir}/download_template_srr.sh
  if [[ $test2 = true ]]; then
   . ${scriptdir}/download_template_srr_parallel.sh
    if [[ $download = true ]]; then
       . ${scriptdir}/download_template_srr_download.sh
    else #download
      echo "******Set download=true when ready to download the fastq files******"; echo ""   
    fi #download      
  if [[ $test3 = true ]]; then
    if [[ $process = true ]]; then
      . ${scriptdir}/download_template_process.sh
    else
      echo "Set process=true to do preprocessing"  
    fi #process
    if [[ $align = true ]]; then
      . ${scriptdir}/download_template_align.sh
    else  
      echo "Set align=true to align files"    
    fi #align
    if [[ $duplicate = true ]]; then
      . ${scriptdir}/download_template_duplicate.sh
    else  
      echo "Set duplicate=true to deduplicate files"     
    fi #duplicate
  else 
    echo "******If the fastq files have been checked with FASTQC, set test3=true" 
  fi #test3
  else 
    echo "******If the SRR FILE is properly formated, set test2 = true to download SRR files******"
  fi #test2
else
  echo "******If the SAMPLE FILE is properly formated, set test1 = true******"
fi #test1

if [[ $test1 = true ]] && [[ $test2 = true ]] && [[ $test3 = true ]] && [[ $test4 = true ]]; then
  if [[ $combine_bam = true ]]; then
    . ${scriptdir}/download_template_bam_combine.sh
  else
    echo "To combine bam files technical replicate set combine_bam=true"
  fi #combine_bam
  if [[ $combine_fastq = true ]]; then
    . ${scriptdir}/download_template_fastq_combine.sh
  else
    echo "To combine fastq files technical replicate set combine_fastq=true"
  fi #combine_fastq
else
  echo "Not combining technical replicates"
fi

if [[ $test1 = true ]] && [[ $test2 = true ]] && [[ $test3 = true ]] && [[ $test4 = true ]] && [[ $test5 = true ]]; then
  if [[ $bio_replicate = true ]]; then
  echo "#--------------Combining biological replicates----------------"
    if [[ $replicate_copy = true ]]; then
      . ${scriptdir}/download_template_replicate_copy.sh    
    else
      echo "To copy files to the proper place to combine biological replicates, set replicate_copy=true"
    fi
    if [[ $replicate_combine = true ]]; then
      . ${scriptdir}/download_template_replicate_copy.sh  
    else
      echo "To combine bam files for biological replicates, set replicate_combine=true"    
    fi       
    . ${scriptdir}/download_template_bio_replicate.sh  
  else
    echo ""; echo "#--------------Not combining biological replicates------------------"
  fi #bio_replicate
  sed '' $infodir/sample_files_${des}.txt; echo "" 
  if [[ $homer = true ]]; then
    echo "HOMER"
  else
    echo "Set homer=true to create tag directories"
  fi
  if [[ $macs2 = true ]]; then
    echo "MACS2"
  else
    echo "Set macs2=true to do peak calling with MACS 2"
  fi      
fi

if [[ $test1 = true ]] && [[ $test2 = true ]] && [[ $test3 = true ]] && [[ $test4 = true ]] && [[ $test5 = true ]] && [[ $macs_14 = true ]]; then
. ${scriptdir}/download_template_macs_14.sh  
##---------------------CLEAN UP MACS14 PARALLEL COMMANDS-----------------------------------------------
#Here are sample sed commands:
#Replace text with NewText
###sed -i -e 's/text/NewText/g' $paralleldir/macs_parallel_${des}.txt
#Delete lines with pattern
###sed -i -e '/pattern/d' $paralleldir/macs_parallel_${des}.txt
#Print lines matching pattern
###sed -i -n '/pattern;/p' $paralleldir/macs_parallel_${des}.txt
sed ' ' $paralleldir/macs_parallel_${des}.txt > $paralleldir/macs_parallel_atac_${des}.txt
sed ' ' $paralleldir/macs_parallel_${des}.txt > $paralleldir/macs_parallel_chip_${des}.txt
sed ' ' $paralleldir/macs_parallel_${des}.txt > $paralleldir/macs_parallel_input_${des}.txt
#sed -i -n '/pattern;/p' $paralleldir/macs_parallel_${des}.txt

echo ""
echo "#------------------CHECK MACS14 PARALLEL COMMANDS------------------------"
sed '' $paralleldir/macs_parallel_${des}.txt
echo "#--------------------------------------------------------"
echo ""
if [[ $macs_14_atac = true ]]; then
  echo "-----------------CHECK MACS14 ATAC Seq PARALLEL COMMANDS-----------------"
  sed '' $paralleldir/macs_parallel_atac_${des}.txt
  echo "#--------------------------------------------------------"; echo ""
fi
if [[ $macs_14_chip = true ]]; then
  echo "-----------------CHECK MACS14 ChIP Seq PARALLEL COMMANDS-----------------"
  sed '' $paralleldir/macs_parallel_chip_${des}.txt
  echo "#--------------------------------------------------------"; echo ""
  echo "-----------------CHECK MACS14 ChIP Seq Input PARALLEL COMMANDS-----------------"
  sed '' $paralleldir/macs_parallel_input_${des}.txt
  echo "#--------------------------------------------------------"; echo ""
fi
echo "Check to make sure MACS14 commands are correct before setting the next option"
echo "Format should be PATH_TO_BAM OUTPUT_NAME fragLength Input"
echo "Transcription factor ChIP-seq should have the full path to the Input file added to the end of each line (separated with a space)"
echo "Set macs_14_chip = true if sample files are transcription factor ChIP-seq with an Input file"
echo "Set macs_14_atac = true if sample files do not have an Input file"
fi #macs_14
