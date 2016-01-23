# The purpose of this script is to look for lineage-specific ATAC-seq peaks.
# First a union is taken of all ATAC-seq peaks for the available cell types
# The union of the peaks is then used to get the original tag coverage (and normalize them) at each location for each cell type ({cell}_Union_Split_Tag_normalize.bedgraph)
# These files are then read into R and a data.frame is constructed
# Cells are split into two groups, an "open" group and a "closed" group
# The peaks that are 1.5 FC greater in the "open" group over the "closed" group are found, if there are more than 100 the list is returned
# This is repeated for all possible combinations of cells in the "open" vs. "closed" group
# 

#Set working directory
wd <- commandArgs(TRUE)
setwd(wd)
#setwd("/mnt/data1/John/Pioneer_Factors")

#Set packages
library(dplyr)
library(pheatmap)
library(RColorBrewer)

#Set options
options(max.print=1000000) 
options(width=900)

#Read in Data
#Set path to bedgraph files
bedgraphPath = "Analysis/ATAC_seq/Union/"
#List the available files
bedgraphFiles <- list.files(path=bedgraphPath, pattern="*_Union_Split_Tag_normalize.bedgraph")
#Loop through the files, read in the data, and assign variables for each file corresponding to the cell type
for (i in bedgraphFiles ) {
  filename <- gsub("_Union_Split_Tag_normalize.bedgraph", "", i)
  cell <- gsub("\\_.*","", i)
  file <- paste0(bedgraphPath,i)
  reads <- read.delim(file, sep = "\t", header=FALSE)
  colnames(reads) <- c("CHR", "Start", "End", cell)
  assign(paste0("reads.",cell), reads)
}

#Create a combined data.frame of all cell types
reads.mat <- subset(reads, select = c(CHR, Start, End))
reads.mat <- cbind(reads.mat, reads.B$B, reads.CD4$CD4, reads.CD8$CD8, reads.CMP$CMP, reads.GMP$GMP, reads.GN$GN, reads.Lsk$Lsk, reads.MEP$MEP, reads.Mono$Mono, reads.NK$NK)
colnames(reads.mat) <- c("CHR", "Start", "End", "B", "CD4", "CD8", "CMP", "GMP", "GN", "Lsk", "MEP", "Mono", "NK")

##Create a matrix of only peaks, add 1, and take the log2, and write this to a file
X <- subset(reads.mat, select = -c(CHR, Start, End))
X.1 <- log2(X+1)
reads.X.1 <- cbind(subset(reads.mat, select = c(CHR, Start, End)), X.1)
write.table(reads.X.1,file="Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#------Test Fold Change 1.5----
fc=1.5
##Unique CD4 Peaks compared to ALL other cell types
CD4.unique.all <- reads.X.1 %>% filter(CD4>0) %>% filter((CD4-B)>fc) %>% filter((CD4-CD8)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc) %>% filter((CD4-NK)>fc) %>% select(CHR, Start, End, CD4)
#------------Take reads.X.1-----Filter CD4>0------------Filter B--------------Filter CD8-----------------Filter CMP------------------Filter GMP----------Filter GN-----------------Filter Lsk----------------Filter MEP--------------Filter Mono---------------Filter NK
write.table(CD4.unique.all,file="Analysis/ATAC_seq/Union/R_output/CD4_unique_all.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

##Unique CD8 Peaks compared to ALL other cell types
CD8.unique.all <- reads.X.1 %>% filter(CD8>0) %>% filter((CD8-B)>fc) %>% filter((CD8-CD4)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((CD8-NK)>fc) %>% select(CHR, Start, End, CD8)
#------------Take reads.X.1-----Filter CD8>0------------Filter B--------------Filter CD4----------------Filter CMP----------------Filter GMP-----------Filter GN-----------------Filter Lsk----------------Filter MEP--------------Filter Mono---------------Filter NK
write.table(CD8.unique.all,file="Analysis/ATAC_seq/Union/R_output/CD8_unique_all.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

##Unique B Peaks compared to ALL other cell types
B.unique.all <- reads.X.1 %>% filter(B>0) %>% filter((B-CD4)>fc) %>% filter((B-CD8)>fc) %>% filter((B-CMP)>fc) %>% filter((B-GMP)>fc) %>% filter((B-GN)>fc) %>% filter((B-Lsk)>fc) %>% filter((B-MEP)>fc) %>% filter((B-Mono)>fc) %>% filter((B-NK)>fc) %>% select(CHR, Start, End, B)
#------------Take reads.X.1-----Filter B>0------------Filter CD4-----------Filter CD8-------------Filter CMP--------------Filter GMP--------Filter GN---------------Filter Lsk--------------Filter MEP-------------Filter Mono---------------Filter NK
write.table(B.unique.all,file="Analysis/ATAC_seq/Union/R_output/B_unique_all.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

##Unique NK Peaks compared to ALL other cell types
NK.unique.all <- reads.X.1 %>% filter(NK>0) %>% filter((NK-B)>fc) %>% filter((NK-CD4)>fc) %>% filter((NK-CD8)>fc) %>% filter((NK-CMP)>fc) %>% filter((NK-GMP)>fc) %>% filter((NK-GN)>fc) %>% filter((NK-Lsk)>fc) %>% filter((NK-MEP)>fc) %>% filter((NK-Mono)>fc) %>% select(CHR, Start, End, NK)
#------------Take reads.X.1-----Filter NK>0----------Filter B-------------Filter CD4----------------Filter CD8--------------Filter CMP----------------Filter GMP----------Filter GN-----------------Filter Lsk----------------Filter MEP--------------Filter Mono---------------Filter NK
write.table(NK.unique.all,file="Analysis/ATAC_seq/Union/R_output/NK_unique_all.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#------Test Fold Change 1----
fc=1
##Unique CD4 Peaks compared to ALL other cell types
CD4.unique.all <- reads.X.1 %>% filter(CD4>0) %>% filter((CD4-B)>fc) %>% filter((CD4-CD8)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc) %>% filter((CD4-NK)>fc) %>% select(CHR, Start, End, CD4)
#------------Take reads.X.1-----Filter CD4>0------------Filter B--------------Filter CD8-----------------Filter CMP------------------Filter GMP----------Filter GN-----------------Filter Lsk----------------Filter MEP--------------Filter Mono---------------Filter NK
write.table(CD4.unique.all,file="Analysis/ATAC_seq/Union/R_output_FC1/CD4_unique_all_FC1.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)
write.table(CD4.unique.all,file="Analysis/ATAC_seq/Union/R_output_Intersect/CD4_unique_all_FC1.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

##Unique CD8 Peaks compared to ALL other cell types
CD8.unique.all <- reads.X.1 %>% filter(CD8>0) %>% filter((CD8-B)>fc) %>% filter((CD8-CD4)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((CD8-NK)>fc) %>% select(CHR, Start, End, CD8)
#------------Take reads.X.1-----Filter CD8>0------------Filter B--------------Filter CD4----------------Filter CMP----------------Filter GMP-----------Filter GN-----------------Filter Lsk----------------Filter MEP--------------Filter Mono---------------Filter NK
write.table(CD8.unique.all,file="Analysis/ATAC_seq/Union/R_output_FC1/CD8_unique_all_FC1.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)
write.table(CD8.unique.all,file="Analysis/ATAC_seq/Union/R_output_Intersect/CD8_unique_all_FC1.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

##Unique B Peaks compared to ALL other cell types
B.unique.all <- reads.X.1 %>% filter(B>0) %>% filter((B-CD4)>fc) %>% filter((B-CD8)>fc) %>% filter((B-CMP)>fc) %>% filter((B-GMP)>fc) %>% filter((B-GN)>fc) %>% filter((B-Lsk)>fc) %>% filter((B-MEP)>fc) %>% filter((B-Mono)>fc) %>% filter((B-NK)>fc) %>% select(CHR, Start, End, B)
#------------Take reads.X.1-----Filter B>0------------Filter CD4-----------Filter CD8-------------Filter CMP--------------Filter GMP--------Filter GN---------------Filter Lsk--------------Filter MEP-------------Filter Mono---------------Filter NK
write.table(B.unique.all,file="Analysis/ATAC_seq/Union/R_output_FC1/B_unique_all_FC1.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)
write.table(B.unique.all,file="Analysis/ATAC_seq/Union/R_output_Intersect/B_unique_all_FC1.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

##Unique NK Peaks compared to ALL other cell types
NK.unique.all <- reads.X.1 %>% filter(NK>0) %>% filter((NK-B)>fc) %>% filter((NK-CD4)>fc) %>% filter((NK-CD8)>fc) %>% filter((NK-CMP)>fc) %>% filter((NK-GMP)>fc) %>% filter((NK-GN)>fc) %>% filter((NK-Lsk)>fc) %>% filter((NK-MEP)>fc) %>% filter((NK-Mono)>fc) %>% select(CHR, Start, End, NK)
#------------Take reads.X.1-----Filter NK>0----------Filter B-------------Filter CD4----------------Filter CD8--------------Filter CMP----------------Filter GMP----------Filter GN-----------------Filter Lsk----------------Filter MEP--------------Filter Mono---------------Filter NK
write.table(NK.unique.all,file="Analysis/ATAC_seq/Union/R_output_FC1/NK_unique_all_FC1.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)
write.table(NK.unique.all,file="Analysis/ATAC_seq/Union/R_output_Intersect/NK_unique_all_FC1.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#------Test Fold Change 0.75--------
fc=0.75
##Unique CD4 Peaks compared to ALL other cell types
CD4.unique.all <- reads.X.1 %>% filter(CD4>0) %>% filter((CD4-B)>fc) %>% filter((CD4-CD8)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc) %>% filter((CD4-NK)>fc) %>% select(CHR, Start, End, CD4)
#------------Take reads.X.1-----Filter CD4>0------------Filter B--------------Filter CD8-----------------Filter CMP------------------Filter GMP----------Filter GN-----------------Filter Lsk----------------Filter MEP--------------Filter Mono---------------Filter NK
write.table(CD4.unique.all,file="Analysis/ATAC_seq/Union/R_output_FC75/CD4_unique_all_FC75.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

##Unique CD8 Peaks compared to ALL other cell types
CD8.unique.all <- reads.X.1 %>% filter(CD8>0) %>% filter((CD8-B)>fc) %>% filter((CD8-CD4)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((CD8-NK)>fc) %>% select(CHR, Start, End, CD8)
#------------Take reads.X.1-----Filter CD8>0------------Filter B--------------Filter CD4----------------Filter CMP----------------Filter GMP-----------Filter GN-----------------Filter Lsk----------------Filter MEP--------------Filter Mono---------------Filter NK
write.table(CD8.unique.all,file="Analysis/ATAC_seq/Union/R_output_FC75/CD8_unique_all_FC75.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#------Test Different Combinations with Fold Change 1.5--------
fc=1.5
#-----Common to NK, CD4, and CD8, but unique to all others---------------
NK.CD4.CD8 <- reads.X.1 %>% 
  filter((CD4-B)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc)  %>% 
  filter((CD8-B)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% 
  filter((NK-B)>fc) %>% filter((NK-CMP)>fc) %>% filter((NK-GMP)>fc) %>% filter((NK-GN)>fc) %>% filter((NK-Lsk)>fc) %>% filter((NK-MEP)>fc) %>% filter((NK-Mono)>fc) %>% select(CHR, Start, End)
write.table(NK.CD4.CD8,file="Analysis/ATAC_seq/Union/R_output_Common/NK_CD4_CD8/NK_CD4_CD8.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#-----Common to NK and CD4, but unique to all others---------------
NK.CD4 <- reads.X.1 %>% 
  filter((CD4-B)>fc) %>% filter((CD4-CD8)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc)  %>% 
  filter((NK-B)>fc) %>% filter((NK-B)>fc) %>% filter((NK-CMP)>fc) %>% filter((NK-GMP)>fc) %>% filter((NK-GN)>fc) %>% filter((NK-Lsk)>fc) %>% filter((NK-MEP)>fc) %>% filter((NK-Mono)>fc) %>% select(CHR, Start, End)
write.table(NK.CD4,file="Analysis/ATAC_seq/Union/R_output_Common/NK_CD4/NK_CD4.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#-----Common to NK and CD8, but unique to all others---------------
NK.CD8 <- reads.X.1 %>% 
  filter((CD8-B)>fc) %>% filter((CD8-CD4)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((NK-B)>fc) %>% 
  filter((NK-CD4)>fc) %>% filter((NK-CMP)>fc) %>% filter((NK-GMP)>fc) %>% filter((NK-GN)>fc) %>% filter((NK-Lsk)>fc) %>% filter((NK-MEP)>fc) %>% filter((NK-Mono)>fc) %>% select(CHR, Start, End)
write.table(NK.CD8,file="Analysis/ATAC_seq/Union/R_output_Common/NK_CD8/NK_CD8.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#-----Common to CD4 and CD8, but unique to all others---------------
#CD4.CD8 <- reads.X.1 %>% filter((CD8-B)>fc) %>% filter((CD8-NK)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((CD4-B)>fc) %>% filter((CD4-NK)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc)  %>% filter((NK-B)<fc) %>% filter((NK-CMP)<fc) %>% filter((NK-GMP)<fc) %>% filter((NK-GN)<fc) %>% filter((NK-Lsk)<fc) %>% filter((NK-MEP)<fc) %>% filter((NK-Mono)<fc) %>% select(CHR, Start, End)
CD4.CD8 <- reads.X.1 %>% 
  filter((CD4-B)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc)  %>% 
  filter((CD8-B)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% 
  filter((NK-B)<fc) %>% filter((NK-CMP)<fc) %>% filter((NK-GMP)<fc) %>% filter((NK-GN)<fc) %>% filter((NK-Lsk)<fc) %>% filter((NK-MEP)<fc) %>% filter((NK-Mono)<fc) %>% select(CHR, Start, End)
write.table(CD4.CD8,file="Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8/CD4_CD8.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#-----Common to CD4 and CD8, but unique to all others, including NK cells---------------
#CD4.CD8 <- reads.X.1 %>% filter((CD8-B)>fc) %>% filter((CD8-NK)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((CD4-B)>fc) %>% filter((CD4-NK)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc) %>% select(CHR, Start, End)
fc=1.5
CD4.CD8 <- reads.X.1 %>% 
  filter((CD4-B)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc)  %>% filter((CD4-NK)>fc) %>%
  filter((CD8-B)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((CD8-NK)>fc) %>% select(CHR, Start, End)
write.table(CD4.CD8,file="Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8_FC_NK/CD4_CD8.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

combinationTesting <- function(fold.change, fold.change.2, cutoff){
  #----Making heatmap----
  cells <- c("B", "CD4", "CD8", "CMP", "GMP", "GN", "Lsk", "MEP", "Mono", "NK")
  positive = NULL
  fc = fold.change
  fc2 = fold.change.2
  testnum = 0
  #----Tests combination of 1-----
  tests <- combn(cells,1)
  ncol(tests)
  for(i in 1:ncol(tests)) 
  {
    open <- tests[,i]  
    closed <- setdiff(cells, open)
    var1.1 <- paste0(open[1],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var1.2 <- paste0(open[1],"-",closed[2],">",fc," & ",open[1],"-",closed[2],"<",fc2)
    var1.3 <- paste0(open[1],"-",closed[3],">",fc," & ",open[1],"-",closed[3],"<",fc2)
    var1.4 <- paste0(open[1],"-",closed[4],">",fc," & ",open[1],"-",closed[4],"<",fc2)
    var1.5 <- paste0(open[1],"-",closed[5],">",fc," & ",open[1],"-",closed[5],"<",fc2)
    var1.6 <- paste0(open[1],"-",closed[6],">",fc," & ",open[1],"-",closed[6],"<",fc2)
    var1.7 <- paste0(open[1],"-",closed[7],">",fc," & ",open[1],"-",closed[7],"<",fc2)
    var1.8 <- paste0(open[1],"-",closed[8],">",fc," & ",open[1],"-",closed[8],"<",fc2)
    var1.9 <- paste0(open[1],"-",closed[9],">",fc," & ",open[1],"-",closed[9],"<",fc2)
    testnum = testnum + 1
    results <- reads.X.1 %>% 
      filter_(var1.1) %>% filter_(var1.2) %>% filter_(var1.3) %>% filter_(var1.4) %>% filter_(var1.5) %>% filter_(var1.6) %>% filter_(var1.7) %>% filter_(var1.8) %>% filter_(var1.9) %>%
      select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
    if (nrow(results) > cutoff) {
      positive = rbind(positive,paste(open[1],sep="."))
      assign(paste(open[1], sep="."), results, envir = .GlobalEnv)
      cell_sys=paste(open[1], sep="_")
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/",cell_sys)
      system(paste0("mkdir -p ",dirFC))
      write.table(results, file=paste0(dirFC,"/",cell_sys,".bed"),sep="\t", col.names = F, row.names = F, quote=FALSE)
    }
  }
  
  #----Tests combination of 2----
  tests <- combn(cells,2)
  ncol(tests)
  for(i in 1:ncol(tests)) 
  {
    open <- tests[,i]  
    closed <- setdiff(cells, open)
    var1.1 <- paste0(open[1],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var1.2 <- paste0(open[1],"-",closed[2],">",fc," & ",open[1],"-",closed[2],"<",fc2)
    var1.3 <- paste0(open[1],"-",closed[3],">",fc," & ",open[1],"-",closed[3],"<",fc2)
    var1.4 <- paste0(open[1],"-",closed[4],">",fc," & ",open[1],"-",closed[4],"<",fc2)
    var1.5 <- paste0(open[1],"-",closed[5],">",fc," & ",open[1],"-",closed[5],"<",fc2)
    var1.6 <- paste0(open[1],"-",closed[6],">",fc," & ",open[1],"-",closed[6],"<",fc2)
    var1.7 <- paste0(open[1],"-",closed[7],">",fc," & ",open[1],"-",closed[7],"<",fc2)
    var1.8 <- paste0(open[1],"-",closed[8],">",fc," & ",open[1],"-",closed[8],"<",fc2)
    var2.1 <- paste0(open[2],"-",closed[1],">",fc," & ",open[2],"-",closed[1],"<",fc2)
    var2.2 <- paste0(open[2],"-",closed[2],">",fc," & ",open[2],"-",closed[2],"<",fc2)
    var2.3 <- paste0(open[2],"-",closed[3],">",fc," & ",open[2],"-",closed[3],"<",fc2)
    var2.4 <- paste0(open[2],"-",closed[4],">",fc," & ",open[2],"-",closed[4],"<",fc2)
    var2.5 <- paste0(open[2],"-",closed[5],">",fc," & ",open[2],"-",closed[5],"<",fc2)
    var2.6 <- paste0(open[2],"-",closed[6],">",fc," & ",open[2],"-",closed[6],"<",fc2)
    var2.7 <- paste0(open[2],"-",closed[7],">",fc," & ",open[2],"-",closed[7],"<",fc2)
    var2.8 <- paste0(open[2],"-",closed[8],">",fc," & ",open[2],"-",closed[8],"<",fc2)
    testnum = testnum + 1
    results <- reads.X.1 %>% 
      filter_(var1.1) %>% filter_(var1.2) %>% filter_(var1.3) %>% filter_(var1.4) %>% filter_(var1.5) %>% filter_(var1.6) %>% filter_(var1.7) %>% filter_(var1.8) %>%
      filter_(var2.1) %>% filter_(var2.2) %>% filter_(var2.3) %>% filter_(var2.4) %>% filter_(var2.5) %>% filter_(var2.6) %>% filter_(var2.7) %>% filter_(var2.8) %>%
      select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
    if (nrow(results) > cutoff) {
      positive = rbind(positive,paste(open[1], open[2],sep="."))
      assign(paste(open[1], open[2], sep="."), results, envir = .GlobalEnv)
      cell_sys=paste(open[1], open[2], sep="_")
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/",cell_sys)
      system(paste0("mkdir -p ",dirFC))
      write.table(results, file=paste0(dirFC,"/",cell_sys,".bed"),sep="\t", col.names = F, row.names = F, quote=FALSE)
    }
  }
  testnum
  
  #----Tests combination of 3----
  tests <- combn(cells,3)
  ncol(tests)
  for(i in 1:ncol(tests)) 
  {
    open <- tests[,i]  
    closed <- setdiff(cells, open)
    var1.1 <- paste0(open[1],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var1.2 <- paste0(open[1],"-",closed[2],">",fc," & ",open[1],"-",closed[2],"<",fc2)
    var1.3 <- paste0(open[1],"-",closed[3],">",fc," & ",open[1],"-",closed[3],"<",fc2)
    var1.4 <- paste0(open[1],"-",closed[4],">",fc," & ",open[1],"-",closed[4],"<",fc2)
    var1.5 <- paste0(open[1],"-",closed[5],">",fc," & ",open[1],"-",closed[5],"<",fc2)
    var1.6 <- paste0(open[1],"-",closed[6],">",fc," & ",open[1],"-",closed[6],"<",fc2)
    var1.7 <- paste0(open[1],"-",closed[7],">",fc," & ",open[1],"-",closed[7],"<",fc2)
    var2.1 <- paste0(open[2],"-",closed[1],">",fc," & ",open[2],"-",closed[1],"<",fc2)
    var2.2 <- paste0(open[2],"-",closed[2],">",fc," & ",open[2],"-",closed[2],"<",fc2)
    var2.3 <- paste0(open[2],"-",closed[3],">",fc," & ",open[2],"-",closed[3],"<",fc2)
    var2.4 <- paste0(open[2],"-",closed[4],">",fc," & ",open[2],"-",closed[4],"<",fc2)
    var2.5 <- paste0(open[2],"-",closed[5],">",fc," & ",open[2],"-",closed[5],"<",fc2)
    var2.6 <- paste0(open[2],"-",closed[6],">",fc," & ",open[2],"-",closed[6],"<",fc2)
    var2.7 <- paste0(open[2],"-",closed[7],">",fc," & ",open[2],"-",closed[7],"<",fc2)
    var3.1 <- paste0(open[3],"-",closed[1],">",fc," & ",open[3],"-",closed[1],"<",fc2)
    var3.2 <- paste0(open[3],"-",closed[2],">",fc," & ",open[3],"-",closed[2],"<",fc2)
    var3.3 <- paste0(open[3],"-",closed[3],">",fc," & ",open[3],"-",closed[3],"<",fc2)
    var3.4 <- paste0(open[3],"-",closed[4],">",fc," & ",open[3],"-",closed[4],"<",fc2)
    var3.5 <- paste0(open[3],"-",closed[5],">",fc," & ",open[3],"-",closed[5],"<",fc2)
    var3.6 <- paste0(open[3],"-",closed[6],">",fc," & ",open[3],"-",closed[6],"<",fc2)
    var3.7 <- paste0(open[3],"-",closed[7],">",fc," & ",open[3],"-",closed[7],"<",fc2)
    testnum = testnum + 1
    results <- reads.X.1 %>% 
      filter_(var1.1) %>% filter_(var1.2) %>% filter_(var1.3) %>% filter_(var1.4) %>% filter_(var1.5) %>% filter_(var1.6) %>% filter_(var1.7) %>% 
      filter_(var2.1) %>% filter_(var2.2) %>% filter_(var2.3) %>% filter_(var2.4) %>% filter_(var2.5) %>% filter_(var2.6) %>% filter_(var2.7) %>% 
      filter_(var3.1) %>% filter_(var3.2) %>% filter_(var3.3) %>% filter_(var3.4) %>% filter_(var3.5) %>% filter_(var3.6) %>% filter_(var3.7) %>% 
      select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
    if (nrow(results) > cutoff) {
      positive = rbind(positive,paste(open[1], open[2], open[3],sep="."))
      assign(paste(open[1], open[2], open[3], sep="."), results, envir = .GlobalEnv)
      cell_sys=paste(open[1], open[2], open[3], sep="_")
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/",cell_sys)
      system(paste0("mkdir -p ",dirFC))
      write.table(results, file=paste0(dirFC,"/",cell_sys,".bed"),sep="\t", col.names = F, row.names = F, quote=FALSE)
    }
  }
  testnum
  
  #----Tests combination of 4----
  tests <- combn(cells,4)
  ncol(tests)
  for(i in 1:ncol(tests)) 
  {
    open <- tests[,i]  
    closed <- setdiff(cells, open)
    var1.1 <- paste0(open[1],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var1.2 <- paste0(open[1],"-",closed[2],">",fc," & ",open[1],"-",closed[2],"<",fc2)
    var1.3 <- paste0(open[1],"-",closed[3],">",fc," & ",open[1],"-",closed[3],"<",fc2)
    var1.4 <- paste0(open[1],"-",closed[4],">",fc," & ",open[1],"-",closed[4],"<",fc2)
    var1.5 <- paste0(open[1],"-",closed[5],">",fc," & ",open[1],"-",closed[5],"<",fc2)
    var1.6 <- paste0(open[1],"-",closed[6],">",fc," & ",open[1],"-",closed[6],"<",fc2)
    var2.1 <- paste0(open[2],"-",closed[1],">",fc," & ",open[2],"-",closed[1],"<",fc2)
    var2.2 <- paste0(open[2],"-",closed[2],">",fc," & ",open[2],"-",closed[2],"<",fc2)
    var2.3 <- paste0(open[2],"-",closed[3],">",fc," & ",open[2],"-",closed[3],"<",fc2)
    var2.4 <- paste0(open[2],"-",closed[4],">",fc," & ",open[2],"-",closed[4],"<",fc2)
    var2.5 <- paste0(open[2],"-",closed[5],">",fc," & ",open[2],"-",closed[5],"<",fc2)
    var2.6 <- paste0(open[2],"-",closed[6],">",fc," & ",open[2],"-",closed[6],"<",fc2)
    var3.1 <- paste0(open[3],"-",closed[1],">",fc," & ",open[3],"-",closed[1],"<",fc2)
    var3.2 <- paste0(open[3],"-",closed[2],">",fc," & ",open[3],"-",closed[2],"<",fc2)
    var3.3 <- paste0(open[3],"-",closed[3],">",fc," & ",open[3],"-",closed[3],"<",fc2)
    var3.4 <- paste0(open[3],"-",closed[4],">",fc," & ",open[3],"-",closed[4],"<",fc2)
    var3.5 <- paste0(open[3],"-",closed[5],">",fc," & ",open[3],"-",closed[5],"<",fc2)
    var3.6 <- paste0(open[3],"-",closed[6],">",fc," & ",open[3],"-",closed[6],"<",fc2)
    var4.1 <- paste0(open[4],"-",closed[1],">",fc," & ",open[4],"-",closed[1],"<",fc2)
    var4.2 <- paste0(open[4],"-",closed[2],">",fc," & ",open[4],"-",closed[2],"<",fc2)
    var4.3 <- paste0(open[4],"-",closed[3],">",fc," & ",open[4],"-",closed[3],"<",fc2)
    var4.4 <- paste0(open[4],"-",closed[4],">",fc," & ",open[4],"-",closed[4],"<",fc2)
    var4.5 <- paste0(open[4],"-",closed[5],">",fc," & ",open[4],"-",closed[5],"<",fc2)
    var4.6 <- paste0(open[4],"-",closed[6],">",fc," & ",open[4],"-",closed[6],"<",fc2)
    testnum = testnum + 1
    results <- reads.X.1 %>% 
      filter_(var1.1) %>% filter_(var1.2) %>% filter_(var1.3) %>% filter_(var1.4) %>% filter_(var1.5) %>% filter_(var1.6) %>% 
      filter_(var2.1) %>% filter_(var2.2) %>% filter_(var2.3) %>% filter_(var2.4) %>% filter_(var2.5) %>% filter_(var2.6) %>% 
      filter_(var3.1) %>% filter_(var3.2) %>% filter_(var3.3) %>% filter_(var3.4) %>% filter_(var3.5) %>% filter_(var3.6) %>% 
      filter_(var4.1) %>% filter_(var4.2) %>% filter_(var4.3) %>% filter_(var4.4) %>% filter_(var4.5) %>% filter_(var4.6) %>% 
      select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
    if (nrow(results) > cutoff) {
      positive = rbind(positive,paste(open[1], open[2], open[3], open[4],sep="."))
      assign(paste(open[1], open[2], open[3], open[4], sep="."), results, envir = .GlobalEnv)
      cell_sys=paste(open[1], open[2], open[3], open[4], sep="_")
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/",cell_sys)
      system(paste0("mkdir -p ",dirFC))
      write.table(results, file=paste0(dirFC,"/",cell_sys,".bed"),sep="\t", col.names = F, row.names = F, quote=FALSE)
    }
  }
  
  #----Tests combination of 5----
  tests <- combn(cells,5)
  ncol(tests)
  for(i in 1:ncol(tests)) 
  {
    open <- tests[,i]  
    closed <- setdiff(cells, open)
    var1.1 <- paste0(open[1],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var1.2 <- paste0(open[1],"-",closed[2],">",fc," & ",open[1],"-",closed[2],"<",fc2)
    var1.3 <- paste0(open[1],"-",closed[3],">",fc," & ",open[1],"-",closed[3],"<",fc2)
    var1.4 <- paste0(open[1],"-",closed[4],">",fc," & ",open[1],"-",closed[4],"<",fc2)
    var1.5 <- paste0(open[1],"-",closed[5],">",fc," & ",open[1],"-",closed[5],"<",fc2)
    var2.1 <- paste0(open[2],"-",closed[1],">",fc," & ",open[2],"-",closed[1],"<",fc2)
    var2.2 <- paste0(open[2],"-",closed[2],">",fc," & ",open[2],"-",closed[2],"<",fc2)
    var2.3 <- paste0(open[2],"-",closed[3],">",fc," & ",open[2],"-",closed[3],"<",fc2)
    var2.4 <- paste0(open[2],"-",closed[4],">",fc," & ",open[2],"-",closed[4],"<",fc2)
    var2.5 <- paste0(open[2],"-",closed[5],">",fc," & ",open[2],"-",closed[5],"<",fc2)
    var3.1 <- paste0(open[3],"-",closed[1],">",fc," & ",open[3],"-",closed[1],"<",fc2)
    var3.2 <- paste0(open[3],"-",closed[2],">",fc," & ",open[3],"-",closed[2],"<",fc2)
    var3.3 <- paste0(open[3],"-",closed[3],">",fc," & ",open[3],"-",closed[3],"<",fc2)
    var3.4 <- paste0(open[3],"-",closed[4],">",fc," & ",open[3],"-",closed[4],"<",fc2)
    var3.5 <- paste0(open[3],"-",closed[5],">",fc," & ",open[3],"-",closed[5],"<",fc2)
    var4.1 <- paste0(open[4],"-",closed[1],">",fc," & ",open[4],"-",closed[1],"<",fc2)
    var4.2 <- paste0(open[4],"-",closed[2],">",fc," & ",open[4],"-",closed[2],"<",fc2)
    var4.3 <- paste0(open[4],"-",closed[3],">",fc," & ",open[4],"-",closed[3],"<",fc2)
    var4.4 <- paste0(open[4],"-",closed[4],">",fc," & ",open[4],"-",closed[4],"<",fc2)
    var4.5 <- paste0(open[4],"-",closed[5],">",fc," & ",open[4],"-",closed[5],"<",fc2)
    var5.1 <- paste0(open[5],"-",closed[1],">",fc," & ",open[5],"-",closed[1],"<",fc2)
    var5.2 <- paste0(open[5],"-",closed[2],">",fc," & ",open[5],"-",closed[2],"<",fc2)
    var5.3 <- paste0(open[5],"-",closed[3],">",fc," & ",open[5],"-",closed[3],"<",fc2)
    var5.4 <- paste0(open[5],"-",closed[4],">",fc," & ",open[5],"-",closed[4],"<",fc2)
    var5.5 <- paste0(open[5],"-",closed[5],">",fc," & ",open[5],"-",closed[5],"<",fc2)
    testnum = testnum + 1
    results <- reads.X.1 %>% 
      filter_(var1.1) %>% filter_(var1.2) %>% filter_(var1.3) %>% filter_(var1.4) %>% filter_(var1.5) %>% 
      filter_(var2.1) %>% filter_(var2.2) %>% filter_(var2.3) %>% filter_(var2.4) %>% filter_(var2.5) %>% 
      filter_(var3.1) %>% filter_(var3.2) %>% filter_(var3.3) %>% filter_(var3.4) %>% filter_(var3.5) %>% 
      filter_(var4.1) %>% filter_(var4.2) %>% filter_(var4.3) %>% filter_(var4.4) %>% filter_(var4.5) %>% 
      filter_(var5.1) %>% filter_(var5.2) %>% filter_(var5.3) %>% filter_(var5.4) %>% filter_(var5.5) %>% 
      select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
    if (nrow(results) > cutoff) {
      positive = rbind(positive,paste(open[1], open[2], open[3], open[4], open[5], sep="."))
      assign(paste(open[1], open[2], open[3], open[4], open[5], sep="."), results, envir = .GlobalEnv)
      cell_sys=paste(open[1], open[2], open[3], open[4], open[5], sep="_")
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/",cell_sys)
      system(paste0("mkdir -p ",dirFC))
      write.table(results, file=paste0(dirFC,"/",cell_sys,".bed"),sep="\t", col.names = F, row.names = F, quote=FALSE)
    }
  }
  
  #----Tests combination of 6----
  tests <- combn(cells,6)
  ncol(tests)
  for(i in 1:ncol(tests)) 
  {
    open <- tests[,i]  
    closed <- setdiff(cells, open)
    var1.1 <- paste0(open[1],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var1.2 <- paste0(open[1],"-",closed[2],">",fc," & ",open[1],"-",closed[2],"<",fc2)
    var1.3 <- paste0(open[1],"-",closed[3],">",fc," & ",open[1],"-",closed[3],"<",fc2)
    var1.4 <- paste0(open[1],"-",closed[4],">",fc," & ",open[1],"-",closed[4],"<",fc2)
    var2.1 <- paste0(open[2],"-",closed[1],">",fc," & ",open[2],"-",closed[1],"<",fc2)
    var2.2 <- paste0(open[2],"-",closed[2],">",fc," & ",open[2],"-",closed[2],"<",fc2)
    var2.3 <- paste0(open[2],"-",closed[3],">",fc," & ",open[2],"-",closed[3],"<",fc2)
    var2.4 <- paste0(open[2],"-",closed[4],">",fc," & ",open[2],"-",closed[4],"<",fc2)
    var3.1 <- paste0(open[3],"-",closed[1],">",fc," & ",open[3],"-",closed[1],"<",fc2)
    var3.2 <- paste0(open[3],"-",closed[2],">",fc," & ",open[3],"-",closed[2],"<",fc2)
    var3.3 <- paste0(open[3],"-",closed[3],">",fc," & ",open[3],"-",closed[3],"<",fc2)
    var3.4 <- paste0(open[3],"-",closed[4],">",fc," & ",open[3],"-",closed[4],"<",fc2)
    var4.1 <- paste0(open[4],"-",closed[1],">",fc," & ",open[4],"-",closed[1],"<",fc2)
    var4.2 <- paste0(open[4],"-",closed[2],">",fc," & ",open[4],"-",closed[2],"<",fc2)
    var4.3 <- paste0(open[4],"-",closed[3],">",fc," & ",open[4],"-",closed[3],"<",fc2)
    var4.4 <- paste0(open[4],"-",closed[4],">",fc," & ",open[4],"-",closed[4],"<",fc2)
    var5.1 <- paste0(open[5],"-",closed[1],">",fc," & ",open[5],"-",closed[1],"<",fc2)
    var5.2 <- paste0(open[5],"-",closed[2],">",fc," & ",open[5],"-",closed[2],"<",fc2)
    var5.3 <- paste0(open[5],"-",closed[3],">",fc," & ",open[5],"-",closed[3],"<",fc2)
    var5.4 <- paste0(open[5],"-",closed[4],">",fc," & ",open[5],"-",closed[4],"<",fc2)
    var6.1 <- paste0(open[6],"-",closed[1],">",fc," & ",open[6],"-",closed[1],"<",fc2)
    var6.2 <- paste0(open[6],"-",closed[2],">",fc," & ",open[6],"-",closed[2],"<",fc2)
    var6.3 <- paste0(open[6],"-",closed[3],">",fc," & ",open[6],"-",closed[3],"<",fc2)
    var6.4 <- paste0(open[6],"-",closed[4],">",fc," & ",open[6],"-",closed[4],"<",fc2)
    testnum = testnum + 1
    results <- reads.X.1 %>% 
      filter_(var1.1) %>% filter_(var1.2) %>% filter_(var1.3) %>% filter_(var1.4) %>% 
      filter_(var2.1) %>% filter_(var2.2) %>% filter_(var2.3) %>% filter_(var2.4) %>% 
      filter_(var3.1) %>% filter_(var3.2) %>% filter_(var3.3) %>% filter_(var3.4) %>% 
      filter_(var4.1) %>% filter_(var4.2) %>% filter_(var4.3) %>% filter_(var4.4) %>% 
      filter_(var5.1) %>% filter_(var5.2) %>% filter_(var5.3) %>% filter_(var5.4) %>% 
      filter_(var6.1) %>% filter_(var6.2) %>% filter_(var6.3) %>% filter_(var6.4) %>% 
      select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
    if (nrow(results) > cutoff) {
      positive = rbind(positive,paste(open[1], open[2], open[3], open[4], open[5], open[6], sep="."))
      assign(paste(open[1], open[2], open[3], open[4], open[5], open[6], sep="."), results, envir = .GlobalEnv)
      cell_sys=paste(open[1], open[2], open[3], open[4], open[5], open[6], sep="_")
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/",cell_sys)
      system(paste0("mkdir -p ",dirFC))
      write.table(results, file=paste0(dirFC,"/",cell_sys,".bed"),sep="\t", col.names = F, row.names = F, quote=FALSE)
    }
  }
  
  #----Tests combination of 7----
  tests <- combn(cells,7)
  ncol(tests)
  for(i in 1:ncol(tests)) 
  {
    open <- tests[,i]  
    closed <- setdiff(cells, open)
    var1.1 <- paste0(open[1],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var1.2 <- paste0(open[1],"-",closed[2],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var1.3 <- paste0(open[1],"-",closed[3],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var2.1 <- paste0(open[2],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var2.2 <- paste0(open[2],"-",closed[2],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var2.3 <- paste0(open[2],"-",closed[3],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var3.1 <- paste0(open[3],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var3.2 <- paste0(open[3],"-",closed[2],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var3.3 <- paste0(open[3],"-",closed[3],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var4.1 <- paste0(open[4],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var4.2 <- paste0(open[4],"-",closed[2],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var4.3 <- paste0(open[4],"-",closed[3],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var5.1 <- paste0(open[5],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var5.2 <- paste0(open[5],"-",closed[2],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var5.3 <- paste0(open[5],"-",closed[3],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var6.1 <- paste0(open[6],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var6.2 <- paste0(open[6],"-",closed[2],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var6.3 <- paste0(open[6],"-",closed[3],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var7.1 <- paste0(open[7],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var7.2 <- paste0(open[7],"-",closed[2],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var7.3 <- paste0(open[7],"-",closed[3],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    testnum = testnum + 1
    results <- reads.X.1 %>% 
      filter_(var1.1) %>% filter_(var1.2) %>% filter_(var1.3) %>% 
      filter_(var2.1) %>% filter_(var2.2) %>% filter_(var2.3) %>% 
      filter_(var3.1) %>% filter_(var3.2) %>% filter_(var3.3) %>% 
      filter_(var4.1) %>% filter_(var4.2) %>% filter_(var4.3) %>% 
      filter_(var5.1) %>% filter_(var5.2) %>% filter_(var5.3) %>% 
      filter_(var6.1) %>% filter_(var6.2) %>% filter_(var6.3) %>% 
      filter_(var7.1) %>% filter_(var7.2) %>% filter_(var7.3) %>% 
      select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
    if (nrow(results) > cutoff) {
      positive = rbind(positive,paste(open[1], open[2], open[3], open[4], open[5], open[6], open[7], sep="."))
      assign(paste(open[1], open[2], open[3], open[4], open[5], open[6], open[7], sep="."), results, envir = .GlobalEnv)
      cell_sys=paste(open[1], open[2], open[3], open[4], open[5], open[6], open[7], sep="_")
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/",cell_sys)
      system(paste0("mkdir -p ",dirFC))
      write.table(results, file=paste0(dirFC,"/",cell_sys,".bed"),sep="\t", col.names = F, row.names = F, quote=FALSE)
    }
  }
  
  #----Tests combination of 8----
  tests <- combn(cells,8)
  ncol(tests)
  for(i in 1:ncol(tests)) 
  {
    open <- tests[,i]  
    closed <- setdiff(cells, open)
    var1.1 <- paste0(open[1],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var1.2 <- paste0(open[1],"-",closed[2],">",fc," & ",open[1],"-",closed[2],"<",fc2)
    var2.1 <- paste0(open[2],"-",closed[1],">",fc," & ",open[2],"-",closed[1],"<",fc2)
    var2.2 <- paste0(open[2],"-",closed[2],">",fc," & ",open[2],"-",closed[2],"<",fc2)
    var3.1 <- paste0(open[3],"-",closed[1],">",fc," & ",open[3],"-",closed[1],"<",fc2)
    var3.2 <- paste0(open[3],"-",closed[2],">",fc," & ",open[3],"-",closed[2],"<",fc2)
    var4.1 <- paste0(open[4],"-",closed[1],">",fc," & ",open[4],"-",closed[1],"<",fc2)
    var4.2 <- paste0(open[4],"-",closed[2],">",fc," & ",open[4],"-",closed[2],"<",fc2)
    var5.1 <- paste0(open[5],"-",closed[1],">",fc," & ",open[5],"-",closed[1],"<",fc2)
    var5.2 <- paste0(open[5],"-",closed[2],">",fc," & ",open[5],"-",closed[2],"<",fc2)
    var6.1 <- paste0(open[6],"-",closed[1],">",fc," & ",open[6],"-",closed[1],"<",fc2)
    var6.2 <- paste0(open[6],"-",closed[2],">",fc," & ",open[6],"-",closed[2],"<",fc2)
    var7.1 <- paste0(open[7],"-",closed[1],">",fc," & ",open[7],"-",closed[1],"<",fc2)
    var7.2 <- paste0(open[7],"-",closed[2],">",fc," & ",open[7],"-",closed[2],"<",fc2)
    var8.1 <- paste0(open[8],"-",closed[1],">",fc," & ",open[8],"-",closed[1],"<",fc2)
    var8.2 <- paste0(open[8],"-",closed[2],">",fc," & ",open[8],"-",closed[2],"<",fc2)
    testnum = testnum + 1
    results <- reads.X.1 %>% 
      filter_(var1.1) %>% filter_(var1.2) %>%
      filter_(var2.1) %>% filter_(var2.2) %>%
      filter_(var3.1) %>% filter_(var3.2) %>%
      filter_(var4.1) %>% filter_(var4.2) %>%
      filter_(var5.1) %>% filter_(var5.2) %>%
      filter_(var6.1) %>% filter_(var6.2) %>%
      filter_(var7.1) %>% filter_(var7.2) %>%
      filter_(var8.1) %>% filter_(var8.2) %>%
      select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
    if (nrow(results) > cutoff) {
      positive = rbind(positive,paste(open[1], open[2], open[3], open[4], open[5], open[6], open[7], open[8], sep="."))
      assign(paste(open[1], open[2], open[3], open[4], open[5], open[6], open[7], open[8], sep="."), results, envir = .GlobalEnv)
      cell_sys=paste(open[1], open[2], open[3], open[4], open[5], open[6], open[7], open[8], sep="_")
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/",cell_sys)
      system(paste0("mkdir -p ",dirFC))
      write.table(results, file=paste0(dirFC,"/",cell_sys,".bed"),sep="\t", col.names = F, row.names = F, quote=FALSE)
    }
  }
  
  #---Tests combination of 9----
  tests <- combn(cells,9)
  ncol(tests)
  for(i in 1:ncol(tests)) 
  {
    open <- tests[,i]  
    closed <- setdiff(cells, open)
    var1.1 <- paste0(open[1],"-",closed[1],">",fc," & ",open[1],"-",closed[1],"<",fc2)
    var2.1 <- paste0(open[2],"-",closed[1],">",fc," & ",open[2],"-",closed[1],"<",fc2)
    var3.1 <- paste0(open[3],"-",closed[1],">",fc," & ",open[3],"-",closed[1],"<",fc2)
    var4.1 <- paste0(open[4],"-",closed[1],">",fc," & ",open[4],"-",closed[1],"<",fc2)
    var5.1 <- paste0(open[5],"-",closed[1],">",fc," & ",open[5],"-",closed[1],"<",fc2)
    var6.1 <- paste0(open[6],"-",closed[1],">",fc," & ",open[6],"-",closed[1],"<",fc2)
    var7.1 <- paste0(open[7],"-",closed[1],">",fc," & ",open[7],"-",closed[1],"<",fc2)
    var8.1 <- paste0(open[8],"-",closed[1],">",fc," & ",open[8],"-",closed[1],"<",fc2)
    var9.1 <- paste0(open[9],"-",closed[1],">",fc," & ",open[9],"-",closed[1],"<",fc2)
    testnum = testnum + 1
    results <- reads.X.1 %>% 
      filter_(var1.1) %>%
      filter_(var2.1) %>%
      filter_(var3.1) %>%
      filter_(var4.1) %>%
      filter_(var5.1) %>%
      filter_(var6.1) %>%
      filter_(var7.1) %>%
      filter_(var8.1) %>%
      filter_(var9.1) %>%
      select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
    if (nrow(results) > cutoff) {
      positive = rbind(positive,paste(open[1], open[2], open[3], open[4], open[5], open[6], open[7], open[8], open[9], sep="."))
      assign(paste(open[1], open[2], open[3], open[4], open[5], open[6], open[7], open[8], open[9], sep="."), results, envir = .GlobalEnv)
      cell_sys=paste(open[1], open[2], open[3], open[4], open[5], open[6], open[7], open[8], open[9], sep="_")
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/",cell_sys)
      system(paste0("mkdir -p ",dirFC))
      write.table(results, file=paste0(dirFC,"/",cell_sys,".bed"),sep="\t", col.names = F, row.names = F, quote=FALSE)
    }
  }
  
  #---Print Results----
  print(testnum)
  for(i in 1:nrow(positive)){
    print(positive[i,])
  }

  return(positive)
}

quantile.range <- quantile(as.matrix(reads.X.1[,4:13]), probs = seq(0, 1, 0.001))
quantile.range

palette.breaks <- seq(quantile.range["5%"], quantile.range["99%"], 0.01)
palette.breaks

#positiveResults <- combinationTesting(0.5, 1.5, 100)
#positiveResults <- combinationTesting(1, 5, 100)
positiveResults <- combinationTesting(1.5, 5, 100)
#positiveResults <- combinationTesting(2, 5, 100)

for(i in 1:nrow(positiveResults)){
  if (i==1){
    labels = nrow(get(positiveResults[1]))
  }
  else{
    labels = c(labels, nrow(get(positiveResults[i])))
  }
}
labels.row1.5 = labels

for(i in 1:nrow(positiveResults)){
  names_D <- positiveResults[i]
  assign(paste0(names_D,".pca"),prcomp(get(names_D)[,4:13], center = TRUE, scale. = TRUE))
  if (i==1){
    heatmap <- get(names_D)
    heatmap.pca <- get(paste0(names_D,".pca"))[[3]]
    rownames <- names_D
  }
  else{
    #---Combine Modules----
    #dplyr::bind_rows(y, z)
    #Append z to y as new rows
    heatmap <- bind_rows(heatmap, get(names_D))
    heatmap.pca <- rbind(heatmap.pca, get(paste0(names_D,".pca"))[[3]])
    rownames <- c(rownames, names_D)
  }
}
heatmap1.5 <- distinct(heatmap)
row.names(heatmap.pca) <- rownames
heatmap1.5.pca <- heatmap.pca

#dplyr::anti_join(a, b, by = "x1")
#All rows in a that do not have a match in b.
#heatmap.remaining <- anti_join(reads.X.1, heatmap1.5)
#heatmap.all <- bind_rows(heatmap, heatmap.remaining)


fc=1.5
bedPath = paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/")
breaks=seq(0, 5, by=0.2)
col.pal  <- colorRampPalette(brewer.pal(9,"Reds"))(length(breaks)-1)
pdf(paste0(bedPath,"Heatmap_All.pdf"), width=7.5, height=11, onefile=FALSE)
pheatmap(heatmap.all[,4:13], cellwidth = 15, scale = "none",
         show_rownames = F, show_colnames = T, breaks=breaks, color=col.pal,
         main = "ATAC-seq Peaks",
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

breaks=seq(0, 5, by=0.2)
col.pal  <- colorRampPalette(brewer.pal(9,"Reds"))(length(breaks)-1)
pdf(paste0(bedPath,"Heatmap_Remaining.pdf"), width=7.5, height=11, onefile=FALSE)
pheatmap(heatmap.remaining[,4:13], cellwidth = 15, scale = "none", kmeans_k = 50,
         show_rownames = F, show_colnames = T, breaks=breaks, color=col.pal,
         main = "ATAC-seq Peaks",
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()


positiveResults <- combinationTesting(0.5, 1.5, 100)

for(i in 1:nrow(positiveResults)){
  if (i==1){
    labels = nrow(get(positiveResults[1]))
  }
  else{
    labels = c(labels, nrow(get(positiveResults[i])))
  }
}
labels.row0.5 = labels
labels.row = c(labels.row1.5,labels.row0.5)

for(i in 1:nrow(positiveResults)){
  names_D <- positiveResults[i]
  assign(paste0(names_D,".pca"),prcomp(get(names_D)[,4:13], center = TRUE, scale. = TRUE))
  if (i==1){
    heatmap <- get(names_D)
    heatmap.pca <- get(paste0(names_D,".pca"))[[3]]
    rownames <- names_D
  }
  else{
    #---Combine Modules----
    #dplyr::bind_rows(y, z)
    #Append z to y as new rows
    heatmap <- bind_rows(heatmap, get(names_D))
    heatmap.pca <- rbind(heatmap.pca, get(paste0(names_D,".pca"))[[3]])
    rownames <- c(rownames, names_D)
  }
}
heatmap0.5 <- distinct(heatmap)
row.names(heatmap.pca) <- rownames
heatmap0.5.pca <- heatmap.pca

#dplyr::anti_join(a, b, by = "x1")
#All rows in a that do not have a match in b.
heatmap.remaining <- anti_join(reads.X.1, heatmap1.5)
heatmap.remaining <- anti_join(heatmap.remaining, heatmap0.5)
#heatmap.all <- bind_rows(heatmap, heatmap.remaining)
heatmap.all <- bind_rows(heatmap1.5, heatmap0.5)
heatmap.pca.all <- bind_rows(data.frame(heatmap1.5.pca), data.frame(heatmap0.5.pca))

fc=1.5

bedPath = paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/")
breaks=seq(0, 5, by=0.2)
col.pal  <- colorRampPalette(brewer.pal(9,"Reds"))(length(breaks)-1)
pdf(paste0(bedPath,"Heatmap_All.pdf"), width=7.5, height=11, onefile=FALSE)
pheatmap(heatmap.all[,4:13], cellwidth = 15, scale = "none",
         show_rownames = F, show_colnames = T, breaks=breaks, color=col.pal,
         main = "ATAC-seq Peaks",
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

pdf(paste0(bedPath,"Heatmap_Reduced.pdf"), width=7.5, height=11, onefile=FALSE)
pheatmap(heatmap.pca.all, cellwidth = 30, scale = "none",
         show_rownames = T, show_colnames = T, 
         main = "ATAC-seq Peaks Reduced",
         cluster_rows = FALSE, cluster_cols = FALSE, labels_row = labels.row, fontsize_row = 8)
dev.off()


#--------------------------------------------
#--------------------------------------------
#--------------------------------------------
#--------------------------------------------
#--------------------------------------------
#--------------------------------------------
for(i in 1:nrow(positiveResults)){
  names_D <- positiveResults[i]
  names_U <- gsub("[.]","_",positiveResults[i])
  bedPath = paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/")
  file <- paste0(bedPath,names_U,"/",names_U,".bed")
  reads <- read.delim(file, sep = "\t", header=FALSE)
  colnames(reads) <- c("CHR", "Start", "End", "B", "CD4", "CD8", "CMP", "GMP", "GN", "Lsk", "MEP", "Mono", "NK")
  
  assign(paste(names_D),reads)
  assign(paste0(names_D,".pca"),prcomp(get(names_D)[,4:13], center = TRUE, scale. = TRUE))
  if (i==1){
    heatmap <- get(names_D)
    heatmap.pca <- get(paste0(names_D,".pca"))[[3]]
    rownames <- names_D
  }
  else{
    #---Combine Modules----
    #dplyr::bind_rows(y, z)
    #Append z to y as new rows
    heatmap <- bind_rows(heatmap, get(names_D))
    heatmap.pca <- rbind(heatmap.pca, get(paste0(names_D,".pca"))[[3]])
    rownames <- c(rownames, names_D)
  }
}
heatmap <- distinct(heatmap)

#dplyr::anti_join(a, b, by = "x1")
#All rows in a that do not have a match in b.
heatmap.remaining <- reads.X.1
rm(heatmap.remaining)
heatmap.remaining <- anti_join(reads.X.1, reads.X.1[1:500000,])
heatmap.all <- bind_rows(heatmap, heatmap.remaining)
row.names(heatmap.pca) <- rownames
colnames(heatmap)

head(CD4.CD8,2)
head(CD4.CD8.old,2)
fc=1.5
CD4.CD8.old <- reads.X.1 %>% 
  filter((CD4-B)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc) %>% filter((CD4-NK)>fc) %>%
  filter((CD8-B)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((CD8-NK)>fc) %>% select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
system("mkdir -p Analysis/ATAC_seq/Union/R_output_Common_2/CD4_CD8")
write.table(CD4.CD8,file="Analysis/ATAC_seq/Union/R_output_Common_2/CD4_CD8/CD4_CD8.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)


#positive <- positiveResults
for(i in 1:nrow(positive)){
  names_D <- positive[i]
  names_U <- gsub("[.]","_",positive[i])
  bedPath = paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/")
  file <- paste0(bedPath,names_U,"/",names_U,".bed")
  reads <- read.delim(file, sep = "\t", header=FALSE)
  colnames(reads) <- c("CHR", "Start", "End", "B", "CD4", "CD8", "CMP", "GMP", "GN", "Lsk", "MEP", "Mono", "NK")
  
  assign(paste(names_D),reads)
  assign(paste0(names_D,".pca"),prcomp(get(names_D)[,4:13], center = TRUE, scale. = TRUE))
  if (i==1){
    heatmap <- get(names_D)
    heatmap.pca <- get(paste0(names_D,".pca"))[[3]]
    rownames <- names_D
  }
  else{
    #---Combine Modules----
    #dplyr::bind_rows(y, z)
    #Append z to y as new rows
    heatmap <- bind_rows(heatmap, get(names_D))
    heatmap.pca <- rbind(heatmap.pca, get(paste0(names_D,".pca"))[[3]])
    rownames <- c(rownames, names_D)
  }
}
heatmap <- distinct(heatmap)

#dplyr::anti_join(a, b, by = "x1")
#All rows in a that do not have a match in b.
heatmap.remaining <- reads.X.1
heatmap.remaining <- anti_join(heatmap.remaining, heatmap)
heatmap.all <- bind_rows(heatmap, heatmap.remaining)
row.names(heatmap.pca) <- rownames
colnames(heatmap)

#----Heatmap---
quantile.range <- quantile(as.matrix(heatmap[,4:13]), probs = seq(0, 1, 0.01))
#palette.breaks <- seq(quantile.range["5%"], quantile.range["99%"], 0.01)
#palette.breaks
#breaks=seq(0, round(tail(quantile.range,2)[1]), by=0.2)

breaks=seq(0, 5, by=0.2)
col.pal  <- colorRampPalette(brewer.pal(9,"Reds"))(length(breaks)-1)
pdf(paste0(bedPath,"Heatmap.pdf"), width=7.5, height=11, onefile=FALSE)
pheatmap(heatmap[,4:13], cellwidth = 15, scale = "none",
         show_rownames = F, show_colnames = T, breaks=breaks, color=col.pal,
         main = "ATAC-seq Peaks",
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

labels.row = nrow(get(positive[1]))
for(i in 2:nrow(positive)){
  labels.row = c(labels.row, nrow(get(positive[i])))
}

pdf(paste0(bedPath,"Heatmap_Reduced.pdf"), onefile=FALSE)
nrow(positiveResults)
pheatmap(heatmap.pca, cellwidth = 30, scale = "none",
         show_rownames = T, show_colnames = T, 
         main = "ATAC-seq Peaks Reduced",
         cluster_rows = TRUE, cluster_cols = FALSE, labels_row = labels.row)
dev.off()


