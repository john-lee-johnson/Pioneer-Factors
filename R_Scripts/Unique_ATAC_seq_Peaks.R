##Set working directory
wd <- commandArgs(TRUE)
setwd(wd)
#setwd("/mnt/data1/John/Pioneer_Factors")
#Read in Data
#datadir="/mnt/data1/VahediLab/PTF_Team/Data"
library(dplyr)
library(pheatmap)
library(WGCNA)
options(max.print=1000000) 
options(width=900)

bedgraphPath = "Analysis/ATAC_seq/Union/"
bedgraphFiles <- list.files(path=bedgraphPath, pattern="*_Union_Split_Tag_normalize.bedgraph")
for (i in bedgraphFiles ) {
  filename <- gsub("_Union_Split_Tag_normalize.bedgraph", "", i)
  cell <- gsub("\\_.*","", i)
  file <- paste0(bedgraphPath,i)
  reads <- read.delim(file, sep = "\t", header=FALSE)
  colnames(reads) <- c("CHR", "Start", "End", cell)
  assign(paste0("reads.",cell), reads)
}
reads.mat <- subset(reads, select = c(CHR, Start, End))
reads.mat <- cbind(reads.mat, reads.B$B, reads.CD4$CD4, reads.CD8$CD8, reads.CMP$CMP, reads.GMP$GMP, reads.GN$GN, reads.Lsk$Lsk, reads.MEP$MEP, reads.Mono$Mono, reads.NK$NK)
colnames(reads.mat) <- c("CHR", "Start", "End", "B", "CD4", "CD8", "CMP", "GMP", "GN", "Lsk", "MEP", "Mono", "NK")

##Create a matrix of only peaks
X <- subset(reads.mat, select = -c(CHR, Start, End))
X.1 <- log2(X+1)
reads.X.1 <- cbind(subset(reads.mat, select = c(CHR, Start, End)), X.1)
write.table(reads.X.1,file="Analysis/ATAC_seq/Union/R_output/B_CD4_CD8_CMP_GMP_GN_Lsk_MEP_Mono_NK_union_log2.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

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

fc=0.75
##Unique CD4 Peaks compared to ALL other cell types
CD4.unique.all <- reads.X.1 %>% filter(CD4>0) %>% filter((CD4-B)>fc) %>% filter((CD4-CD8)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc) %>% filter((CD4-NK)>fc) %>% select(CHR, Start, End, CD4)
#------------Take reads.X.1-----Filter CD4>0------------Filter B--------------Filter CD8-----------------Filter CMP------------------Filter GMP----------Filter GN-----------------Filter Lsk----------------Filter MEP--------------Filter Mono---------------Filter NK
write.table(CD4.unique.all,file="Analysis/ATAC_seq/Union/R_output_FC75/CD4_unique_all_FC75.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

##Unique CD8 Peaks compared to ALL other cell types
CD8.unique.all <- reads.X.1 %>% filter(CD8>0) %>% filter((CD8-B)>fc) %>% filter((CD8-CD4)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((CD8-NK)>fc) %>% select(CHR, Start, End, CD8)
#------------Take reads.X.1-----Filter CD8>0------------Filter B--------------Filter CD4----------------Filter CMP----------------Filter GMP-----------Filter GN-----------------Filter Lsk----------------Filter MEP--------------Filter Mono---------------Filter NK
write.table(CD8.unique.all,file="Analysis/ATAC_seq/Union/R_output_FC75/CD8_unique_all_FC75.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

fc=1.5
#-----Common to NK, CD4, and CD8, but unique to all others---------------
NK.CD4.CD8 <- reads.X.1 %>% 
  filter((CD4-B)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc)  %>% 
  filter((CD8-B)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% 
  filter((NK-B)>fc) %>% filter((NK-CMP)>fc) %>% filter((NK-GMP)>fc) %>% filter((NK-GN)>fc) %>% filter((NK-Lsk)>fc) %>% filter((NK-MEP)>fc) %>% filter((NK-Mono)>fc) %>% select(CHR, Start, End)
write.table(NK.CD4.CD8,file="Analysis/ATAC_seq/Union/R_output_Common/NK_CD4_CD8/NK_CD4_CD8.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#-----Common to NK and CD4, but unique to all others---------------
NK.CD4 <- reads.X.1 %>% 
  filter((CD4-B)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc)  %>% 
  filter((CD8-B)<fc) %>% filter((CD8-CMP)<fc) %>% filter((CD8-GMP)<fc) %>% filter((CD8-GN)<fc) %>% filter((CD8-Lsk)<fc) %>% filter((CD8-MEP)<fc) %>% filter((CD8-Mono)<fc) %>% 
  filter((NK-B)>fc) %>% filter((NK-CMP)>fc) %>% filter((NK-GMP)>fc) %>% filter((NK-GN)>fc) %>% filter((NK-Lsk)>fc) %>% filter((NK-MEP)>fc) %>% filter((NK-Mono)>fc) %>% select(CHR, Start, End)
write.table(NK.CD4,file="Analysis/ATAC_seq/Union/R_output_Common/NK_CD4/NK_CD4.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#-----Common to NK and CD8, but unique to all others---------------
NK.CD8 <- reads.X.1 %>% 
  filter((CD4-B)<fc) %>% filter((CD4-CMP)<fc) %>% filter((CD4-GMP)<fc) %>% filter((CD4-GN)<fc) %>% filter((CD4-Lsk)<fc) %>% filter((CD4-MEP)<fc) %>% filter((CD4-Mono)<fc)  %>% 
  filter((CD8-B)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((NK-B)>fc) %>% 
  filter((NK-CMP)>fc) %>% filter((NK-GMP)>fc) %>% filter((NK-GN)>fc) %>% filter((NK-Lsk)>fc) %>% filter((NK-MEP)>fc) %>% filter((NK-Mono)>fc) %>% select(CHR, Start, End)
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
CD4.CD8 <- reads.X.1 %>% 
  filter((CD4-B)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc)  %>% filter((CD4-NK)>fc) %>%
  filter((CD8-B)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((CD8-NK)>fc) %>% select(CHR, Start, End)
write.table(CD4.CD8,file="Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8_FC_NK/CD4_CD8.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)


row1 <- c(15, 2, 2, 2, 2, 2, 2, 2, 2, 2)
row2 <- c(2, 17, 2, 2, 2, 2, 2, 2, 2, 2)
row3 <- c(2, 2, 15, 2, 2, 2, 2, 2, 2, 2)
row4 <- c(2, 2, 2, 21, 2, 2, 2, 2, 2, 2)
row5 <- c(2, 2, 2, 2, 16, 2, 2, 2, 2, 2)
row6 <- c(2, 2, 2, 2, 2, 15, 2, 2, 2, 2)
row7 <- c(2, 2, 2, 2, 2, 2, 16, 2, 2, 2)
row8 <- c(2, 2, 2, 2, 2, 2, 2, 17, 2, 2)
row9 <- c(2, 2, 2, 2, 2, 2, 2, 2, 18, 2)
row10 <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 19)
row11 <- c(2, 16, 2, 2, 2, 2, 2, 2, 2, 19)
row12 <- c(2, 2, 15, 2, 2, 2, 2, 2, 2, 19)
row13 <- c(2, 17, 19, 2, 3, 2, 2, 2, 2, 19)
mat <- rbind(row1, row2, row3, row4, row5, row6, row7, row8, row9, row10, row11, row12, row13)
colnames(mat) <- c("B","CD4","CD8","CMP","GMP","GN","Lsk","MEP","Mono","NK")
rownames(mat) <- c("B", "CD4", "CD8", "CMP", "GMP", "GN", "Lsk", "MEP", "Mono", "NK", "CD4.NK1", "CD8.NK1", "CD4.CD8.NK")
mat
#Lsk.B.Mono <- reads.X.1 %>% filter((Lsk-CD4)>fc) %>% filter((Lsk-CD8)>fc) %>% filter((Lsk-NK)>fc) %>% filter((B-CD4)>fc) %>% filter((B-CD8)>fc) %>% filter((B-NK)>fc) %>% filter((Mono-CD4)>fc) %>% filter((Mono-NK)>fc) %>% filter((Mono-CD8)>fc) %>% select(CHR, Start, End) #%>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((CD4-B)>fc) %>% filter((CD4-NK)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc) %>% select(CHR, Start, End)
#write.table(Lsk.B.Mono,file="Analysis/ATAC_seq/Union/R_output_Close/NK_CD4_CD8/Lsk_B_Mono_close.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)
#Lsk.B.Mono <- reads.X.1 %>% filter((Lsk-CD4)>fc) %>% filter((Lsk-CD8)>fc) %>% filter((Lsk-NK)>fc) %>% filter((Lsk-B)>fc) %>% filter((Mono-CD4)>fc) %>% filter((Mono-NK)>fc) %>% filter((Mono-CD8)>fc) %>% select(CHR, Start, End) #%>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((CD4-B)>fc) %>% filter((CD4-NK)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc) %>% select(CHR, Start, End)
#write.table(Lsk.B.Mono,file="Analysis/ATAC_seq/Union/R_output_Close/NK_CD4_CD8/Lsk_B_Mono_close.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

?pheatmap


# A 2-dimensional dataset
df <- rbind(matrix(rnorm(50, sd = 0.3), ncol = 5), matrix(rnorm(50, mean = 1.5, sd = 0.3), ncol = 5), matrix(rnorm(50, mean = 5.5, sd = 0.3), ncol = 5))
colnames(df) <- c("B", "CD4", "CD8", "NK", "Lsk")
# k-means
(cl <- kmeans(reads.X.1[,4:13], centers=25, iter.max = 10, nstart = 25))
plot(df, col = cl$cluster)
points(cl$centers, col = 1:5, pch = 8)
# determine the best number of clusters
y=numeric(25);
for(i in 1:25) {
  cl=kmeans(df, i)
  y[i]=100*(1-cl$tot.withinss/cl$totss);
}
plot(y, xlab="Number of clusters", ylab="% of explained variance", type='b')


