##Set working directory
wd <- commandArgs(TRUE)
setwd(wd)
#getwd()
#setwd("/mnt/data1/John/Pioneer_Factors")
#setwd("~/Desktop/John/Pioneer_Factors")
setwd("/mnt/data1/John/Pioneer-Factors")
#Read in Data
#datadir="/mnt/data1/VahediLab/PTF_Team/Data"
library(rafalib)
library(genefilter)
library(gplots)
library(dplyr)
library(NMF)
library(hexbin)
library(ggplot2)
library(RColorBrewer)
library(circular)
options(max.print=1000000) 
options(width=900)

reads <- read.table("Analysis/ATAC_seq/B_CD4_CD8_Lsk_NK_bam_count_normalized.bedgraph")
colnames(reads)=c("CHR", "Start", "End", "B", "CD4", "CD8", "Lsk", "NK")

##Create a matrix of only peaks
X <- subset(reads, select = -c(CHR, Start, End))
X.1 <- log2(X+1)
reads.X.1 <- cbind(subset(reads, select = c(CHR, Start, End)), X.1)
write.table(reads.X.1,file="Analysis/ATAC_seq/R_output/B_CD4_CD8_Lsk_NK_bam_log2.bedgraph",sep="\t", col.names = F, row.names = F, quote=FALSE)

fc=1.5
##Unique CD4 Peaks with Fold Change
CD4.unique.Lsk_FC <- reads.X.1 %>% filter(CD4>0) %>% filter((CD4-Lsk)>fc) %>% select(CHR, Start, End, CD4)
write.table(CD4.unique.Lsk_FC,file="Analysis/ATAC_seq/R_output/CD4_unique_Lsk.bedgraph",sep="\t", col.names = F, row.names = F, quote=FALSE)
CD4.unique.Lsk_FC <- reads.X.1 %>% filter(CD4>0) %>% filter((CD4-Lsk)>fc) %>% select(CHR, Start, End, B, CD4, CD8, Lsk, NK)
##Unique CD4 Peaks with Fold Change
CD8.unique.Lsk_FC <- reads.X.1 %>% filter(CD8>0) %>% filter((CD8-Lsk)>fc) %>% select(CHR, Start, End, CD8)
write.table(CD8.unique.Lsk_FC,file="Analysis/ATAC_seq/R_output/CD8_unique_Lsk.bedgraph",sep="\t", col.names = F, row.names = F, quote=FALSE)
CD8.unique.Lsk_FC <- reads.X.1 %>% filter(CD8>0) %>% filter((CD8-Lsk)>fc) %>% select(CHR, Start, End, B, CD4, CD8, Lsk, NK)
##Unique NK Peaks with Fold Change
NK.unique.Lsk_FC <- reads.X.1 %>% filter(NK>0) %>% filter((NK-Lsk)>fc) %>% select(CHR, Start, End, NK)
write.table(NK.unique.Lsk_FC,file="Analysis/ATAC_seq/R_output/NK_unique_Lsk_FC.bedgraph",sep="\t", col.names = F, row.names = F, quote=FALSE)
NK.unique.Lsk_FC <- reads.X.1 %>% filter(NK>0) %>% filter((NK-Lsk)>fc) %>% select(CHR, Start, End, B, CD4, CD8, Lsk, NK)
##Unique B Peaks with Fold Change
B.unique.Lsk_FC <- reads.X.1 %>% filter(B>0) %>% filter((B-Lsk)>fc) %>% select(CHR, Start, End, B)
write.table(B.unique.Lsk_FC,file="Analysis/ATAC_seq/R_output/B_unique_Lsk_FC.bedgraph",sep="\t", col.names = F, row.names = F, quote=FALSE)
B.unique.Lsk_FC <- reads.X.1 %>% filter(B>0) %>% filter((B-Lsk)>fc) %>% select(CHR, Start, End, B, CD4, CD8, Lsk, NK)

recombine <- merge(B.unique.Lsk_FC, CD4.unique.Lsk_FC, all=TRUE)
recombine <- merge(recombine, CD8.unique.Lsk_FC, all=TRUE)
recombine <- merge(recombine, NK.unique.Lsk_FC, all=TRUE)

##Unique CD4 Peaks compared to All other cell types
CD4.unique.all <- recombine %>% filter(CD4>0) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-CD8)>fc) %>% filter((CD4-B)>fc)  %>% filter((CD4-NK)>fc)  %>% select(CHR, Start, End, CD4)
CD4.unique.all <- recombine %>% filter(CD4>0) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-CD8)>fc) %>% filter((CD4-B)>fc)  %>% filter((CD4-NK)>fc)  %>% select(CHR, Start, End, B, CD4, CD8, NK, Lsk)
write.table(CD4.unique.all,file="Analysis/ATAC_seq/R_output/CD4_unique_all.bedgraph",sep="\t", col.names = F, row.names = F, quote=FALSE)
#CD4.unique.all <- recombine %>% filter(CD4>0) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-CD8)>fc) %>% filter((CD4-B)>fc)  %>% filter((CD4-NK)>fc)  %>% select(CHR, Start, End, B, CD4, CD8, Lsk, NK)
##Unique CD8 Peaks compared to All other cell types
CD8.unique.all <- recombine %>% filter(CD8>0) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-CD4)>fc) %>% filter((CD8-B)>fc)  %>% filter((CD8-NK)>fc)  %>% select(CHR, Start, End, CD8)
CD8.unique.all <- recombine %>% filter(CD8>0) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-CD4)>fc) %>% filter((CD8-B)>fc)  %>% filter((CD8-NK)>fc)  %>% select(CHR, Start, End, B, CD4, CD8, NK, Lsk)
write.table(CD8.unique.all,file="Analysis/ATAC_seq/R_output/CD8_unique_all.bedgraph",sep="\t", col.names = F, row.names = F, quote=FALSE)
#CD8.unique.all <- recombine %>% filter(CD8>0) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-CD4)>fc) %>% filter((CD8-B)>fc)  %>% filter((CD8-NK)>fc)  %>% select(CHR, Start, End, B, CD4, CD8, Lsk, NK)
##Unique NK Peaks compared to All other cell types
NK.unique.all <- recombine %>% filter(NK>0) %>% filter((NK-Lsk)>fc) %>% filter((NK-CD8)>fc) %>% filter((NK-CD4)>fc)  %>% filter((NK-CD8)>fc)  %>% select(CHR, Start, End, NK)
NK.unique.all <- recombine %>% filter(NK>0) %>% filter((NK-Lsk)>fc) %>% filter((NK-CD8)>fc) %>% filter((NK-CD4)>fc)  %>% filter((NK-CD8)>fc)  %>% select(CHR, Start, End, B, CD4, CD8, NK, Lsk)
write.table(NK.unique.all,file="Analysis/ATAC_seq/R_output/NK_unique_all.bedgraph",sep="\t", col.names = F, row.names = F, quote=FALSE)
#NK.unique.all <- recombine %>% filter(NK>0) %>% filter((NK-Lsk)>fc) %>% filter((NK-CD8)>fc) %>% filter((NK-CD4)>fc)  %>% filter((NK-CD8)>fc)  %>% select(CHR, Start, End, B, CD4, CD8, Lsk, NK)
##Unique B Peaks compared to All other cell types
B.unique.all <- recombine %>% filter(B>0) %>% filter((B-Lsk)>fc) %>% filter((B-CD4)>fc) %>% filter((B-CD8)>fc)  %>% filter((B-NK)>fc) %>% select(CHR, Start, End, B)
B.unique.all <- recombine %>% filter(B>0) %>% filter((B-Lsk)>fc) %>% filter((B-CD4)>fc) %>% filter((B-CD8)>fc)  %>% filter((B-NK)>fc) %>% select(CHR, Start, End, B, CD4, CD8, Lsk, NK, Lsk)
write.table(B.unique.all,file="Analysis/ATAC_seq/R_output/B_unique_all.bedgraph",sep="\t", col.names = F, row.names = F, quote=FALSE)
#B.unique.all <- recombine %>% filter(B>0) %>% filter((B-Lsk)>fc) %>% filter((B-CD4)>fc) %>% filter((B-CD8)>fc)  %>% filter((B-NK)>fc) %>% select(CHR, Start, End, B, CD4, CD8, Lsk, NK)

recombine.all <- merge(B.unique.all, CD4.unique.all, all=TRUE)
recombine.all <- merge(recombine.all, CD8.unique.all, all=TRUE)
recombine.all <- merge(recombine.all, NK.unique.all, all=TRUE)
recombine.all.minus <- subset(recombine.all, select =-c(CHR, Start, End))
write.table(recombine.all,file="Analysis/ATAC_seq/R_output/Combine_all_unique_peaks.bedgraph",sep="\t", col.names = F, row.names = F, quote=FALSE)
recombine.all <- NULL
recombine.all <- rbind(recombine.all, NK.unique.all)
recombine.all <- rbind(recombine.all, CD4.unique.all)
recombine.all <- rbind(recombine.all, CD8.unique.all)
recombine.all <- rbind(recombine.all, B.unique.all)
recombine.all.minus <- subset(recombine.all, select =-c(CHR, Start, End))


#Heatmap
hmcol <- colorRampPalette(brewer.pal(9, "Purples"))(50)
#png(filename="Analysis/ATAC_seq/R_output/ATAC_seq_Unique_Peaks_Heatmap.png", width=2000, height=1556, units="px", res=300)
#mypar(1,1)
#aheatmap(recombine.all.minus, color = hmcol, scale = "row", Rowv=FALSE, Colv=NA, legend=TRUE, main=paste("Unique ATAC Seq Peaks + nRow = ", nrow(recombine.all.minus)))
postscript(file="Figures/ATAC_seq_Unique_Peaks_Heatmap.eps", onefile=FALSE, horizontal=FALSE)
#aheatmap(recombine.all.minus, color = hmcol, scale = "row", cellwidth=50, Rowv=NA, Colv=NA, legend=TRUE)
aheatmap(recombine.all.minus, color = hmcol, scale = "row", cellwidth=50, Rowv=NA, Colv=NA, legend=TRUE, main=paste("Unique ATAC Seq Peaks + nRow = ", nrow(recombine.all.minus)))
dev.off()
#dev.off()

#Common NK, CD4, and CD8 Peaks with Fold Change > Lsk and B
CD4.CD8.NK.common.B.Lsk_FC <- reads.X.1 %>% filter(CD4>0) %>% filter(CD8>0) %>% filter(NK>0) %>% filter((CD4-Lsk)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((NK-Lsk)>fc) %>% filter((CD4-B)>fc) %>% filter((CD8-B)>fc) %>% filter((NK-B)>fc) %>% select(CHR, Start, End, CD4, CD8, NK, B, Lsk)
write.table(CD4.CD8.NK.common.B.Lsk_FC,file="Analysis/ATAC_seq/R_output/CD4_CD8_NK_common_B_Lsk_FC.bedgraph",sep="\t", col.names = F, row.names = F, quote=FALSE)
CD4.CD8.NK.common.B.Lsk_FC.minus <- subset(CD4.CD8.NK.common.B.Lsk_FC, select =-c(CHR, Start, End))
postscript(file="Figures/CD4_CD8_NK_common_Heatmap.eps", onefile=FALSE, horizontal=FALSE)
png(filename="Figures/CD4_CD8_NK_common_Heatmap.png",width=2000, height=1556, units="px", res=300)
aheatmap(CD4.CD8.NK.common.B.Lsk_FC.minus, color = hmcol, scale = "row", cellwidth=50, Rowv=FALSE, Colv=NA, legend=TRUE, main=paste("Common Accessible CD4, CD8, NK over Lsk, B + nRow = ", nrow(CD4.CD8.NK.common.B.Lsk_FC)))
dev.off()

