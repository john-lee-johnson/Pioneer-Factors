##Set working directory
wd <- commandArgs(TRUE)
setwd(wd)
#setwd("/mnt/data1/John/Pioneer_Factors")

#----Load libraries
library(dplyr)
library(pheatmap)
library(RColorBrewer)

#---Set Options
options(max.print=1000000) 
options(width=900)

#----Normalize coverage files
#Read in file that contains total read numbers per sample
reads.counts <- read.delim("../Pioneer-Factors/info/atac_read_count.txt", sep = " ", header=FALSE)
colnames(reads.counts) <- c("Cell", "Seq", "Investigator", "Reads")
##Set the path
bedgraphPath = "Analysis/ATAC_seq/Union/R_output_merge_coverage/"
##Get the files
bedgraphFiles <- list.files(path=bedgraphPath, pattern="*_union_merge_split_sorted_coverage.bedgraph")
for (i in bedgraphFiles ) {
  ##Get the filename (remove quotes)
  filename <- gsub("_union_merge_split_sorted_coverage.bedgraph", "", i)
  cell <- gsub("\\_.*","", i)
  file <- paste0(bedgraphPath,i)
  reads <- read.delim(file, sep = "\t", header=FALSE)
  reads$V4 <- reads$V4/(reads.counts[reads.counts$Cell==cell,4])
  colnames(reads) <- c("CHR", "Start", "End", cell)
  assign(paste0("reads.",cell), reads)
  write.table(get(paste0("reads.",cell)),file=paste0(bedgraphPath,cell,"_union_merge_split_sorted_coverage_normalized.bedgraph"),sep="\t", col.names = F, row.names = F, quote=FALSE)
}

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
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_Merge",fc,"/",cell_sys)
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
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_Merge",fc,"/",cell_sys)
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
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_Merge",fc,"/",cell_sys)
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
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_Merge",fc,"/",cell_sys)
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
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_Merge",fc,"/",cell_sys)
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
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_Merge",fc,"/",cell_sys)
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
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_Merge",fc,"/",cell_sys)
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
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_Merge",fc,"/",cell_sys)
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
      dirFC <- paste0("Analysis/ATAC_seq/Union/R_output_Common_Merge",fc,"/",cell_sys)
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
rm(fc)

positiveResults <- combinationTesting(0.5, 1.5, 100)
positiveResults <- combinationTesting(1.5, 5, 100)
positiveResults <- combinationTesting(2, 5, 100)