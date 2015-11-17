##Set working directory
wd <- commandArgs(TRUE)
setwd(wd)
library(rtracklayer)
options(max.print=1000000) 
options(width=900)
reads <- read.delim("Analysis/ATAC_seq/B_CD4_CD8_Lsk_NK_merged_common.bed", sep = "\t", header=FALSE)
colnames(reads) <- c("CHR", "Start", "End", "B", "CD4", "CD8", "Lsk", "NK")

write.table(reads[,c(1,2,3,4)],file="Analysis/ATAC_seq/B_union_split.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)
write.table(reads[,c(1,2,3,5)],file="Analysis/ATAC_seq/CD4_union_split.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)
write.table(reads[,c(1,2,3,6)],file="Analysis/ATAC_seq/CD8_union_split.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)
write.table(reads[,c(1,2,3,7)],file="Analysis/ATAC_seq/Lsk_union_split.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)
write.table(reads[,c(1,2,3,8)],file="Analysis/ATAC_seq/NK_union_split.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)