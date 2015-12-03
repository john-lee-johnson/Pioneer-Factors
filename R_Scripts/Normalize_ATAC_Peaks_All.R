##Set working directory
wd <- commandArgs(TRUE)
setwd(wd)
library(rtracklayer)
options(max.print=1000000) 
options(width=900)
#setwd("/mnt/data1/John/Pioneer_Factors")
reads.counts <- read.delim("../Pioneer-Factors/info/atac_read_count.txt", sep = " ", header=FALSE)
colnames(reads.counts) <- c("Cell", "Seq", "Investigator", "Reads")
bedgraphPath = "Analysis/ATAC_seq/Union/"
bedgraphFiles <- list.files(path=bedgraphPath, pattern="*_Union_Split_Tag.bedgraph")
for (i in bedgraphFiles ) {
  filename <- gsub("_Union_Split_Tag.bedgraph", "", i)
  cell <- gsub("\\_.*","", i)
  file <- paste0(bedgraphPath,i)
  reads <- read.delim(file, sep = "\t", header=FALSE)
  reads$V4 <- reads$V4/(reads.counts[reads.counts$Cell==cell,4])
  colnames(reads) <- c("CHR", "Start", "End", cell)
  assign(paste0("reads.",cell), reads)
  write.table(get(paste0("reads.",cell)),file=paste0(bedgraphPath,filename,"_Union_Split_Tag_normalize.bedgraph"),sep="\t", col.names = F, row.names = F, quote=FALSE)
}