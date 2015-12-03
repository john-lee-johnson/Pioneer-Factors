##Set working directory
wd <- commandArgs(TRUE)
setwd(wd)
library(rtracklayer)
options(max.print=1000000) 
options(width=900)

reads.counts <- read.delim("../Pioneer-Factors/info/atac_read_count.txt", sep = " ", header=FALSE)
colnames(reads.counts) <- c("Cell", "Seq", "Investigator", "Reads")

#reads <- read.delim("Analysis/ATAC_seq/B_CD4_CD8_Lsk_NK_bam_count_all.bedgraph", sep = "\t", header=FALSE)
#colnames(reads) <- c("CHR", "Start", "End", "B", "CD4", "CD8", "Lsk", "NK")
#reads$B <- reads$B/(reads.counts[reads.counts$Cell=="B",2])
#reads$CD4 <- reads$CD4/(reads.counts[reads.counts$Cell=="CD4",2])
#reads$CD8 <- reads$CD8/(reads.counts[reads.counts$Cell=="CD8",2])
#reads$Lsk <- reads$Lsk/(reads.counts[reads.counts$Cell=="Lsk",2])
#reads$NK <- reads$NK/(reads.counts[reads.counts$Cell=="NK",2])

removeMax <- function(data,top){
  qnt <- quantile(data, probs=c(.25, top), na.rm = TRUE)
  indx <- arrayInd(which(data>qnt[2]), dim(data))
  return(indx)
}

X <- subset(reads, select = -c(CHR, Start, End))
indx <- removeMax(X,.999)
write.table(reads,file="Analysis/ATAC_seq/B_CD4_CD8_Lsk_NK_bam_count_normalized.bedgraph",sep="\t", col.names = F, row.names = F, quote=FALSE)
