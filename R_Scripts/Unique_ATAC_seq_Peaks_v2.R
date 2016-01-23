# The purpose of this script is to look for lineage-specific ATAC-seq peaks.
# First a union is taken of all ATAC-seq peaks for the available cell types
# The coordinates of the union of the peaks is then used to get the original tag coverage (and normalize them) at each location for each cell type ({cell}_Union_Split_Tag_normalize.bedgraph)
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
library(dendextend)

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

#----Read in Normalized Coverage Files
bedgraphPath = "Analysis/ATAC_seq/Union/R_output_merge_coverage/"
bedgraphFiles <- list.files(path=bedgraphPath, pattern="*_union_merge_split_sorted_coverage_normalized.bedgraph")
reads.X.1 <- read.delim("Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed", sep = "\t", header=FALSE)

for (i in bedgraphFiles ) {
  filename <- gsub("_union_merge_split_sorted_coverage_normalized.bedgraph", "", i)
  cell <- gsub("\\_.*","", i)
  file <- paste0(bedgraphPath,i)
  reads <- read.delim(file, sep = "\t", header=FALSE)
  colnames(reads) <- c("CHR", "Start", "End", cell)
  assign(paste0("merge.",cell), reads)
}
#---Combine normalized coverage into reads.X.1
reads.X.1 <- read.delim("Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2_merge.bed", sep = "\t", header=FALSE)
reads.X.1[,4] <- merge.B[,4]; reads.X.1[,5] <- merge.CD4[,4]; reads.X.1[,6] <- merge.CD8[,4]; reads.X.1[,7] <- merge.CMP[,4]; reads.X.1[,8] <- merge.GMP[,4]; reads.X.1[,9] <- merge.GN[,4]; reads.X.1[,10] <- merge.Lsk[,4]; reads.X.1[,11] <- merge.MEP[,4]; reads.X.1[,12] <- merge.Mono[,4]; reads.X.1[,13] <- merge.NK[,4]
colnames(reads.X.1) <- c("CHR", "Start", "End", "B", "CD4", "CD8", "CMP", "GMP", "GN", "Lsk", "MEP", "Mono", "NK")
#---Add 1 and Take the log2 of the reads
reads.mat <- reads.X.1
X <- subset(reads.mat, select = -c(CHR, Start, End))
X.1 <- log2(X+1)
reads.X.1 <- cbind(subset(reads.mat, select = c(CHR, Start, End)), X.1)

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
  positiveResults <- positive
  return(positive)
}

#Generate the labels for the heatmap
labelsGenerate <- function(positiveResults){
  for(i in 1:nrow(positiveResults)){
    if (i==1){
      labels = nrow(get(positiveResults[1]))
    }
    else{
      labels = c(labels, nrow(get(positiveResults[i])))
    }
  }
  return(labels)
}

#Do PCA to reduce heatmaps
pcaGenerate <- function(positiveResults){
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
  List <- list(heatmap, heatmap.pca, rownames)
  return(List)
}

#Subset rows to highest variance
rowVariance <- function(dat, upper){
  #Removing rows with value > 5
  df <- dat[!apply(dat[,4:13], 1, function(x) {any(x > upper)}),]
  # Find the desired quantile breaks for the entire matrix
  # First get row variances with this function
  rv = apply(as.matrix(df[,4:13]), 1, var)
  # Return the variance cutoff of the 80% percentile
  qt <- quantile(as.matrix(rv) , probs = c(0.2,0.80) )
  #  Next get a logical vector of the rows that have any values outside these breaks
  rows <- apply(as.matrix(rv) , 1 , function(x) ( x > qt[2] ) )
  #  Subset on this vector
  df <- df[ rows , ]
  return(df)
}

#Find optimal number of clusters
clusterPlotting <- function(data){
  #Determine number of clusters
  mydata <- data[,4:13]
  gc()
  wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  for (i in 1:50) wss[i] <- sum(kmeans(mydata,iter.max=10,
                                       centers=i)$withinss)
  plot(1:50, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
}

# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pheatmap(heatmap.all, clustering_callback = callback)

#Generate PDF and return clusters
heatmapClusterReturnFirst <- function(data, fc, maxclust, seedSet){
  set.seed(seedSet)
  basedir = paste0("Analysis/ATAC_seq/Union/R_output_Common_Merge",fc)
  filename <- paste0("Heatmap_Remaining_After_",fc,"_Kmeans_Cluster", maxclust, ".pdf")
  outfile <- paste(basedir,filename,sep="/")
  main <- paste("result for ", maxclust, " clusters", sep="")
  set.seed(seedSet)
  hmx.parameters <- list(data[,4:13],
                         cellwidth = 15,
                         scale = "none",
                         kmeans_k = maxclust,
                         show_rownames = T,
                         show_colnames = T,
                         breaks=breaks,
                         color=col.pal,
                         main = "ATAC-seq Peaks",
                         cluster_rows = TRUE,
                         cluster_cols = FALSE)
  
  # To store cluster mappings and draw
  set.seed(seedSet)
  kmean.hm <- do.call("pheatmap", hmx1.parameters)
  
  # To store cluster mappings and draw
  set.seed(seedSet)
  kmean.hm <- do.call("pheatmap", hmx.parameters)
  
  # To draw on screen 
  set.seed(seedSet)
  #do.call("pheatmap", hmx.parameters)
  
  # To draw to file 
  set.seed(seedSet)
  do.call("pheatmap", c(hmx.parameters, filename=outfile))
  
  # add cluster number to matrix and save
  clustnum <- kmean.hm[["kmeans"]][["cluster"]]
  clustered.data <- cbind(data, clustnum)
  last <- ncol(clustered.data)-1
  
  # Prepare output directory
  basedir = paste0("Analysis/ATAC_seq/Union/R_output_Common_Merge",fc,"/Clusters")
  system(paste0("mkdir -p ",basedir))
  
  for (clust in 1:maxclust){
    filename <- paste("Clusters.for.", clust, ".of.", maxclust, "-clusters.bed", sep="")
    outfile <- paste(basedir,filename,sep="/")
    cluster <- subset(clustered.data, clustered.data$clustnum==clust)[1:last]
    write.table(cluster, row.names=FALSE, col.names=FALSE, file=outfile, quote = FALSE, sep="\t")
  }  
  return(kmean.hm)
}

#Generate PDF and return clusters
heatmapClusterReturn <- function(data, fc, maxclust, seedSet){
  set.seed(seedSet)
  basedir = paste0("Analysis/ATAC_seq/Union/R_output_Common_Merge",fc)
  filename <- paste0("Heatmap_Remaining_After_",fc,"_Kmeans_Cluster", maxclust, ".pdf")
  outfile <- paste(basedir,filename,sep="/")
  main <- paste("result for ", maxclust, " clusters", sep="")
  set.seed(seedSet)
  hmx.parameters <- list(data[,4:13],
                       cellwidth = 15,
                       scale = "none",
                       kmeans_k = maxclust,
                       show_rownames = T,
                       show_colnames = T,
                       breaks=breaks,
                       color=col.pal,
                       main = "ATAC-seq Peaks",
                       cluster_rows = cl2,
                       cluster_cols = FALSE,
                       treeheight_row = 25
                       )
  
  # To store cluster mappings and draw
  set.seed(seedSet)
  kmean.hm <- do.call("pheatmap", hmx1.parameters)
  
  # To store cluster mappings and draw
  set.seed(seedSet)
  kmean.hm <- do.call("pheatmap", hmx.parameters)
  
  # To draw on screen 
  set.seed(seedSet)
  #do.call("pheatmap", hmx.parameters)
  
  # To draw to file 
  set.seed(seedSet)
  do.call("pheatmap", c(hmx.parameters, filename=outfile))
  
  # add cluster number to matrix and save
  clustnum <- kmean.hm[["kmeans"]][["cluster"]]
  clustered.data <- cbind(data, clustnum)
  last <- ncol(clustered.data)-1
  
  # Prepare output directory
  basedir = paste0("Analysis/ATAC_seq/Union/R_output_Common",fc,"/Clusters")
  system(paste0("mkdir -p ",basedir))
  
  for (clust in 1:maxclust){
    filename <- paste("Clusters.for.", clust, ".of.", maxclust, "-clusters.bed", sep="")
    outfile <- paste(basedir,filename,sep="/")
    cluster <- subset(clustered.data, clustered.data$clustnum==clust)[1:last]
    #write.table(cluster, row.names=FALSE, col.names=FALSE, file=outfile, quote = FALSE, sep="\t")
  }  
}

#Examine the data spread
(quantile.range <- quantile(as.matrix(reads.X.1[,4:13]), probs = seq(0, 1, 0.001)))
breaks=seq(0, 4, by=0.2)
breaks=seq(0, 5, by=0.5)
#breaks=seq(0, 6, by=0.5)
length(breaks)
#Set the color palatte
length(breaks)
col.pal  <- colorRampPalette(brewer.pal(8,"Reds"))(length(breaks)-1)
col.pal  <- colorRampPalette(brewer.pal(8,"Reds"))(8)

#Do K-Means Clustering Before Applying Any Subsetting
#Subset on the rows with the highest variance (top 20%)
heatmap.all <- rowVariance(reads.X.1, 5)
clusterPlotting(heatmap.all)
#Generate the PDF and clusters for the first time to observe dendrogram
#Stores as object for later manipulation of row order
gc()
set.seed(124)
kmean.hm1 <- heatmapClusterReturnFirst(heatmap.all,0,15,124)
#Get the hclust object
cl <- kmean.hm1[["tree_row"]]
#Look at the order
cl[["order"]]
#Rotate the dendrogram (The function heatmapClusterReturn will use the order of cl2)
#So make sure the order of cl2 is correct
cl2 <- rotate(cl,as.character(c(3,14,6,12,7,11,1,15,2,13,8,10,4,5,9)))
cl2[["order"]]
set.seed(124)
heatmapClusterReturn(heatmap.all,0,15,124)



heatmapClusterReturn(heatmap.all,0,20,1)

#Apply the first test: a fold change greater than 1.5 and less than 5 (the max)
positiveResults <- combinationTesting(1.5, 5, 100)

#Generate the heatmaps for Fold Change 1.5
labels.row1.5 <- labelsGenerate(positiveResults)
#Get back the rows for the heatmap
listReturn <- pcaGenerate(positiveResults)
#Heatmap for FC 1.5
heatmap1.5 <- distinct(listReturn[[1]])
#Heatmap reduced by PCA for FC 1.5
heatmap1.5.pca <- listReturn[[2]]
#Get the row names for the heatmap
row.names(heatmap1.5.pca) <- listReturn[[3]]

#dplyr::anti_join(a, b, by = "x1")
#All rows in a that do not have a match in b.
#Remove the rows from above
heatmap.remaining <- anti_join(reads.X.1, heatmap1.5)
heatmap.remaining.var  <- rowVariance(heatmap.remaining)
#Generate the PDF and return the coordinates for each cluster
gc()
set.seed(123)
heatmapClusterReturn(heatmap.remaining.var,1.5,15,123)

#Append the remaining rows onto heatmap1,5
#heatmap.all <- bind_rows(heatmap1.5, heatmap.remaining)

#Heatmap of only peaks with greater than 1.5 FC
fc=1.5
bedPath = paste0("Analysis/ATAC_seq/Union/R_output_Common_",fc,"/")
pdf(paste0(bedPath,"Heatmap_1.5_Only.pdf"), width=7.5, height=11, onefile=FALSE)
pheatmap(heatmap1.5[,4:13], cellwidth = 15, scale = "none",
         show_rownames = F, show_colnames = T, breaks=breaks, color=col.pal,
         main = "ATAC-seq Peaks",
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

#Heatmap of only peaks with greater than 1.5 FC reduced by PCA
pdf(paste0(bedPath,"Heatmap_Reduced.pdf"), width=7.5, height=11, onefile=FALSE)
pheatmap(heatmap1.5.pca, cellwidth = 30, scale = "none",
         show_rownames = T, show_colnames = T, 
         main = "ATAC-seq Peaks Reduced",
         cluster_rows = FALSE, cluster_cols = FALSE, labels_row = labels.row1.5, fontsize_row = 8)
dev.off()

#Heatmap of remaining peaks, clustered by k-means
#Subset to the rows with the highest variance (top 20%)
heatmap.remaining <- rowVariance(heatmap.remaining)
#Generate the PDF and return the clusters
heatmapClusterReturn(heatmap.remaining,1.5,15)


#Get the results from the second test: a FC greater than 0.5, but less than 1.5
positiveResults <- combinationTesting(0.5, 1.5, 100)
labels.row1 <- labelsGenerate(positiveResults)
listReturn <- pcaGenerate(positiveResults)
heatmap1 <- distinct(listReturn[[1]])
heatmap1.pca <- listReturn[[2]]
row.names(heatmap1.pca) <- listReturn[[3]]


#dplyr::anti_join(a, b, by = "x1")
#All rows in a that do not have a match in b.
heatmap.remaining <- anti_join(reads.X.1, heatmap1.5)
heatmap.remaining <- anti_join(heatmap.remaining, heatmap1)
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


