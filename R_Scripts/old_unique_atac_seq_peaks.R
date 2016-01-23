
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

reads.mat2 <- subset(reads, select = c(CHR, Start, End))
reads.mat2 <- cbind(reads.mat2, reads.NK$NK, reads.CD4$CD4, reads.CD8$CD8, reads.CMP$CMP, reads.GMP$GMP, reads.GN$GN, reads.Lsk$Lsk, reads.MEP$MEP, reads.Mono$Mono, reads.B$B)
colnames(reads.mat2) <- c("CHR", "Start", "End", "NK", "CD4", "CD8", "CMP", "GMP", "GN", "Lsk", "MEP", "Mono", "B")
X2 <- subset(reads.mat2, select = -c(CHR, Start, End))
X.2 <- log2(X2+1)
reads.X.2 <- cbind(subset(reads.mat2, select = c(CHR, Start, End)), X.2)
colnames(reads.X.2)
head(reads.X.2)
# plot 5-50 kmean clusters
library("pheatmap")
library("RColorBrewer")
col.pal <- brewer.pal(9,"Blues")
# define metrics for clustering
drows1 <- "euclidean"
dcols1 <- "euclidean"
gc()
fc=1.5

#Generate all possible combinations of cell types
cells <- c("B", "CD4", "CD8", "CMP", "GMP", "GN", "Lsk", "MEP", "Mono", "NK")
combn(cells[1:10], 1)
combn(cells[1:10], 2)
combn(cells[1:10], 3)
combn(cells[1:10], 4)
combn(cells[1:10], 5)
combn(cells[1:10], 6)
combn(cells[1:10], 7)
combn(cells[1:10], 8)
combn(cells[1:10], 9)
combn(cells[1:10], 10)

combn(letters[1:4], 2)
combn(letters[1:4], 3)

NK.open <- reads.X.2 %>% filter((NK-Lsk)>fc)
CD4.open <- reads.X.2 %>% filter((CD4-Lsk)>fc)
CD8.open <- reads.X.2 %>% filter((CD8-Lsk)>fc)
CMP.open <- reads.X.2 %>% filter((CMP-Lsk)>fc)
GMP.open <- reads.X.2 %>% filter((GMP-Lsk)>fc)
GN.open <- reads.X.2 %>% filter((GN-Lsk)>fc)
MEP.open <- reads.X.2 %>% filter((MEP-Lsk)>fc)
Mono.open <- reads.X.2 %>% filter((Mono-Lsk)>fc)
B.open <- reads.X.2 %>% filter((B-Lsk)>fc)
merge.reads.X.2 <- merge(NK.open, CD4.open, by = intersect(names(NK.open), names(CD4.open)), all = TRUE)
merge.reads.X.2 <- merge(merge.reads.X.2, CD8.open, by = intersect(names(merge.reads.X.2), names(CD8.open)), all = TRUE)
merge.reads.X.2 <- merge(merge.reads.X.2, CMP.open, by = intersect(names(merge.reads.X.2), names(CMP.open)), all = TRUE)
merge.reads.X.2 <- merge(merge.reads.X.2, GMP.open, by = intersect(names(merge.reads.X.2), names(GMP.open)), all = TRUE)
merge.reads.X.2 <- merge(merge.reads.X.2, GN.open, by = intersect(names(merge.reads.X.2), names(GN.open)), all = TRUE)
merge.reads.X.2 <- merge(merge.reads.X.2, MEP.open, by = intersect(names(merge.reads.X.2), names(MEP.open)), all = TRUE)
merge.reads.X.2 <- merge(merge.reads.X.2, Mono.open, by = intersect(names(merge.reads.X.2), names(Mono.open)), all = TRUE)
merge.reads.X.2 <- merge(merge.reads.X.2, B.open, by = intersect(names(merge.reads.X.2), names(B.open)), all = TRUE)
options(stringsAsFactors = FALSE);

# here we define the adjacency matrix using soft thresholding with beta=6
datExpr <- t(as.matrix(merge.reads.X.2[,4:13]))
gsg <- goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

ADJ1=abs(cor(t(datExpr),use="p"))^6

merge.var <- var(t(merge.reads.X.2[1:5,4:13]))
merge.var
dim(merge.var)


# When you have relatively few genes (<5000) use the following code
k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=datExpr,power=20) 
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")


head(merge.reads.X.2)
head(NK.open)
head(CD4.open)
names(NK.open[,1:3])
head(Lsk.open)
#%>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc) %>% filter((CD4-NK)>fc) %>% select(CHR, Start, End, CD4)

#coordinates <- reads.X.2[reads.X.2$Lsk-reads.X.2$CD4<1.5,]
#coordinates <- coordinates[coordinates$CD4>1.5,]
#coordinates <- coordinates[coordinates$CD8>1.5,1:3]
coordinates <- reads.X.2 %>% filter((Lsk-CD4)>fc) %>% filter((Lsk-CD8)>fc)
head(coordinates)
data<-reads.X.2[reads.X.2$Lsk<3.5,]
data<-data[data$CD4>1.5,]
data<-data[data$CD8>1.5,4:13]
# following code limits the lowest and highest color to 5%, and 95% of your range, respectively
quantile.range <- quantile(as.matrix(data), probs = seq(0, 1, 0.01))
quantile.range
palette.breaks <- seq(quantile.range["5%"], quantile.range["99%"], 0.1)
palette.breaks

#breaks for the core of the distribution
breaks=seq(0, 4, by=0.2) #41 values
breaks=c(breaks,17)
#now add outliers
col.pal  <- colorRampPalette(brewer.pal(9,"Blues"))(length(breaks)-1)
maxclust=100
maxclust=50
for (maxclust in seq(from = 15, to = 50, by = 5)){
  basedir=paste0("Analysis/ATAC_seq/Union/R_output_Discover/Cluster_Plots/K.",maxclust)
  filename <- paste("heatmap.for.", maxclust, "-clusters.pdf", sep="")
  outfile <- paste(basedir,filename,sep="/")
  main <- paste("result for ", maxclust, " clusters", sep="")
  hmx.parameters <- list(data, 
                         color = col.pal,
                         cellwidth = 15, cellheight = 12, scale = "none",
                         kmeans_k = maxclust,
                         breaks=breaks,
                         show_rownames = T, show_colnames = T, 
                         main = main,
                         clustering_method = "average",
                         cluster_rows = FALSE, cluster_cols = FALSE)
  
  # To store cluster mappings and draw
  kmean.hm <- do.call("pheatmap", hmx.parameters)
  
  # To draw to file 
  do.call("pheatmap", c(hmx.parameters, filename=outfile))
  
  # add cluster number to matrix and save
  clustnum <- kmean.hm[["kmeans"]][["cluster"]]
  clustered.data <- cbind(coordinates, clustnum)
  for (clust in 1:maxclust){
    filename <- paste("Clustered.for.", clust, ".of.", maxclust, "-clusters.bed", sep="")
    outfile <- paste(basedir,filename,sep="/")
    write.table(clustered.data[clustered.data$clustnum==clust,1:3], row.names=FALSE, col.names=FALSE, file=outfile, quote = FALSE, sep="\t")
  }  
  gc()
}
#Finished 

MNase <- read.delim("Analysis/MNase_seq/DP/histfile.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
NK.me2.CD4 <- read.delim("Analysis/ChIP_seq/homer/NK/peak_annotation/NK_CD4_H3K4me2_output.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
NK.me3.CD4 <- read.delim("Analysis/ChIP_seq/homer/NK/peak_annotation/NK_CD4_H3K4me3_output.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
NK.27ac.CD4 <- read.delim("Analysis/ChIP_seq/homer/NK/peak_annotation/NK_CD4_H3K27Ac_output.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

dev.off()
y1=0.3
y2=0.6
plot(MNase[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Histone Marks for CD4 at NK Unique Peaks", xlab="Distance From Center", ylab="Tag Coverage",lwd=2.5)
par(new = TRUE)
plot(NK.me2.CD4[,c("Distance","Coverage")], type="l", xlim=range(-2500,2500), ylim=range(c(y1,y2)), col = "blue", xlab="", ylab="",lwd=2.5)
par(new = TRUE)
plot(NK.me3.CD4[,c("Distance","Coverage")], type="l", xlim=range(-2500,2500), ylim=range(c(y1,y2)), col = "orange", xlab="", ylab="",lwd=2.5)
par(new = TRUE)
plot(NK.27ac.CD4[,c("Distance","Coverage")], type="l", xlim=range(-2500,2500), ylim=range(c(y1,y2)), col = "green", xlab="", ylab="",lwd=2.5)
legend("topright", xjust=1, c("H3K4me1","H3K4me2","H3K4me3","H3K27Ac"), lty=c(1,1,1,1), lwd=c(2.5,2.5,2.5,2.5),col=c("red","blue","orange","green"))


#---------------------------------------------------
for(i in 1:nrow(positive)){
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
  }
  else{
    #---Combine Modules----
    #dplyr::bind_rows(y, z)
    #Append z to y as new rows
    heatmap <- bind_rows(heatmap, get(names_D))
    heatmap.pca <- rbind(heatmap.pca, get(paste0(names_D,".pca"))[[3]])
  }
}
heatmap <- distinct(heatmap)

#----Heatmap---
quantile.range <- quantile(as.matrix(heatmap[,4:13]), probs = seq(0, 1, 0.01))
#palette.breaks <- seq(quantile.range["5%"], quantile.range["99%"], 0.01)
#palette.breaks
breaks=seq(0, round(tail(quantile.range,2)[1]), by=0.2)
col.pal  <- colorRampPalette(brewer.pal(9,"Reds"))(length(breaks)-1)
pdf(paste0(bedPath,"Heatmap.pdf"), width=7.5, height=11, onefile=FALSE)
pheatmap(heatmap[,4:13], cellwidth = 15, scale = "none",
         show_rownames = F, show_colnames = T, breaks=breaks, color=col.pal,
         main = "Title",
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

labels.row = nrow(get(positiveResults[1]))
for(i in 2:nrow(positiveResults)){
  labels.row = c(labels.row, nrow(get(positiveResults[i])))
}
labels.row

pdf(paste0(bedPath,"Heatmap_Reduced.pdf"), onefile=FALSE)
pheatmap(heatmap.pca, cellwidth = 30, cellheight = 30, scale = "none",
         show_rownames = T, show_colnames = T, 
         main = "ATAC-seq Peaks Reduced",
         cluster_rows = FALSE, cluster_cols = FALSE, labels_row = labels.row)
dev.off()





#breaks for the core of the distribution
breaks=seq(0, 4, by=0.2)
#now add outliers
breaks=c(breaks,9)
col.pal  <- colorRampPalette(brewer.pal(9,"Reds"))(length(breaks)-1)

pdf(paste0(bedPath,"Heatmap.pdf"), width=7.5, height=11, onefile=FALSE)
pheatmap(heatmap[,4:13], cellwidth = 15, scale = "none",
         show_rownames = F, show_colnames = T, breaks=breaks, color=col.pal,
         main = "ATAC-seq Peak",
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()



#-----PCA on the modules to reduce dimensions-----
CD4.CD8.pca <- prcomp(CD4.CD8[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- CD4.CD8.pca[[3]]
CD4.CD8.NK.pca <- prcomp(CD4.CD8.NK[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, CD4.CD8.NK.pca[[3]])
CD4.pca <- prcomp(CD4[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, CD4.pca[[3]])
CD8.pca <- prcomp(CD8[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, CD8.pca[[3]])
B.CD4.CD8.NK.pca <- prcomp(B.CD4.CD8.NK[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, B.CD4.CD8.NK.pca[[3]])
B.pca <- prcomp(B[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, B.pca[[3]])
NK.pca <- prcomp(NK[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, NK.pca[[3]])
CMP.pca <- prcomp(CMP[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, CMP.pca[[3]])
Lsk.pca <- prcomp(Lsk[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, Lsk.pca[[3]])
Mono.pca <- prcomp(Mono[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, Mono.pca[[3]])
B.Lsk.pca <- prcomp(B.Lsk[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, B.Lsk.pca[[3]])
GN.Mono.pca <- prcomp(GN.Mono[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, GN.Mono.pca[[3]])
B.CD4.CD8.Lsk.NK.pca <- prcomp(B.CD4.CD8.Lsk.NK[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, B.CD4.CD8.Lsk.NK.pca[[3]])
B.CD4.CD8.CMP.GMP.GN.Lsk.Mono.NK.pca <- prcomp(B.CD4.CD8.CMP.GMP.GN.Lsk.Mono.NK[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, B.CD4.CD8.CMP.GMP.GN.Lsk.Mono.NK.pca[[3]])


for(i in 1:nrow(positiveResults)){
  print(positiveResults[i,])
}

positiveResults <- combinationTesting(1.25, 100)

for(i in 1:nrow(positiveResults)){
  print(positiveResults[i,])
}

positiveResults <- combinationTesting(1, 100)

for(i in 1:nrow(positiveResults)){
  print(positiveResults[i,])
}

positiveResults <- combinationTesting(1.25, 100)

for(i in 1:nrow(positiveResults)){
  print(positiveResults[i,])
}


fc=1.5

#---Combine Modules----
#dplyr::bind_rows(y, z)
#Append z to y as new rows
heatmap <- CD4.CD8
heatmap <- bind_rows(heatmap, CD4.CD8.NK)
heatmap <- bind_rows(heatmap, CD4)
heatmap <- bind_rows(heatmap, CD8)
heatmap <- bind_rows(heatmap, B.CD4.CD8.NK)
heatmap <- bind_rows(heatmap, B)
heatmap <- bind_rows(heatmap, NK)
heatmap <- bind_rows(heatmap, CMP)
heatmap <- bind_rows(heatmap, Lsk)
heatmap <- bind_rows(heatmap, Mono)
heatmap <- bind_rows(heatmap, B.Lsk)
heatmap <- bind_rows(heatmap, GN.Mono)
heatmap <- bind_rows(heatmap, B.CD4.CD8.Lsk.NK)
heatmap <- bind_rows(heatmap, B.CD4.CD8.CMP.GMP.GN.Lsk.Mono.NK)
heatmap <- distinct(heatmap)

#dplyr::anti_join(a, b, by = "x1")
#All rows in a that do not have a match in b.
heatmap.remaining <- reads.X.1
heatmap.remaining <- anti_join(heatmap.remaining, heatmap)

nrow(heatmap)-(nrow(reads.X.1)-nrow(heatmap.remaining))

#----Heatmap---
quantile.range <- quantile(as.matrix(heatmap[,4:13]), probs = seq(0, 1, 0.01))
quantile.range
palette.breaks <- seq(quantile.range["5%"], quantile.range["99%"], 0.1)
palette.breaks

#breaks for the core of the distribution
breaks=seq(0, 4, by=0.2)
#now add outliers
breaks=c(breaks,9)
col.pal  <- colorRampPalette(brewer.pal(9,"Reds"))(length(breaks)-1)

pdf('Analysis/ATAC_seq/Union/R_output_Common/Heatmap_Test_1.pdf', width=7.5, height=11, onefile=FALSE)
pheatmap(heatmap[,4:13], cellwidth = 15, scale = "none",
         show_rownames = F, show_colnames = T, breaks=breaks, color=col.pal,
         main = "Title",
         cluster_rows = FALSE, cluster_cols = FALSE)

dev.off()

#-----PCA on the modules to reduce dimensions-----
CD4.CD8.pca <- prcomp(CD4.CD8[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- CD4.CD8.pca[[3]]
CD4.CD8.NK.pca <- prcomp(CD4.CD8.NK[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, CD4.CD8.NK.pca[[3]])
CD4.pca <- prcomp(CD4[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, CD4.pca[[3]])
CD8.pca <- prcomp(CD8[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, CD8.pca[[3]])
B.CD4.CD8.NK.pca <- prcomp(B.CD4.CD8.NK[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, B.CD4.CD8.NK.pca[[3]])
B.pca <- prcomp(B[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, B.pca[[3]])
NK.pca <- prcomp(NK[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, NK.pca[[3]])
CMP.pca <- prcomp(CMP[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, CMP.pca[[3]])
Lsk.pca <- prcomp(Lsk[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, Lsk.pca[[3]])
Mono.pca <- prcomp(Mono[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, Mono.pca[[3]])
B.Lsk.pca <- prcomp(B.Lsk[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, B.Lsk.pca[[3]])
GN.Mono.pca <- prcomp(GN.Mono[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, GN.Mono.pca[[3]])
B.CD4.CD8.Lsk.NK.pca <- prcomp(B.CD4.CD8.Lsk.NK[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, B.CD4.CD8.Lsk.NK.pca[[3]])
B.CD4.CD8.CMP.GMP.GN.Lsk.Mono.NK.pca <- prcomp(B.CD4.CD8.CMP.GMP.GN.Lsk.Mono.NK[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, B.CD4.CD8.CMP.GMP.GN.Lsk.Mono.NK.pca[[3]])

#----Plot Heatmap with Reduced Dimensions----

#---Create Labels for Rows----
labels.row = nrow(get(positiveResults[1]))
for(i in 2:nrow(positiveResults)){
  labels.row = c(labels.row, nrow(get(positiveResults[i])))
}
labels.row

pdf('Analysis/ATAC_seq/Union/R_output_Common/Heatmap_Test_1_Reduced.pdf', onefile=FALSE)
pheatmap(heatmap.pca, cellwidth = 30, cellheight = 30, scale = "none",
         show_rownames = T, show_colnames = T, 
         main = "ATAC-seq Peaks",
         cluster_rows = FALSE, cluster_cols = FALSE, labels_row = labels.row)
dev.off()
?pdf

positiveResults <- combinationTesting(2)

for(i in 1:nrow(positiveResults)){
  print(positiveResults[i,])
}
fc=2
#--------FC2---------------------
#-----FC2 Unique B Peaks compared to ALL other cell types----
B <- reads.X.1 %>% filter((B-CD4)>1.5 & (B-CD4)<3) %>% filter((B-CD8)>fc) %>% filter((B-CMP)>fc) %>% filter((B-GMP)>fc) %>% filter((B-GN)>fc) %>% filter((B-Lsk)>fc) %>% filter((B-MEP)>fc) %>% filter((B-Mono)>fc) %>% filter((B-NK)>fc) %>% select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
system("mkdir -p Analysis/ATAC_seq/Union/R_output_Common_2/B")
write.table(B,file="Analysis/ATAC_seq/Union/R_output_Common_2/B/B.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#----FC2 Unique Lsk Peaks compared to ALL other cell types----
Lsk <- reads.X.1 %>% filter((Lsk-B)>fc) %>% filter((Lsk-CD4)>fc) %>% filter((Lsk-CD8)>fc) %>% filter((Lsk-CMP)>fc) %>% filter((Lsk-GMP)>fc) %>% filter((Lsk-GN)>fc) %>% filter((Lsk-MEP)>fc) %>% filter((Lsk-Mono)>fc) %>% filter((Lsk-NK)>fc) %>% select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
system("mkdir -p Analysis/ATAC_seq/Union/R_output_Common_2/Lsk")
write.table(Lsk,file="Analysis/ATAC_seq/Union/R_output_Common_2/Lsk/Lsk.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#----FC2 Unique NK Peaks compared to ALL other cell types----
NK <- reads.X.1 %>% filter((NK-B)>fc) %>% filter((NK-CD4)>fc) %>% filter((NK-CD8)>fc) %>% filter((NK-CMP)>fc) %>% filter((NK-GMP)>fc) %>% filter((NK-GN)>fc) %>% filter((NK-Lsk)>fc) %>% filter((NK-MEP)>fc) %>% filter((NK-Mono)>fc) %>% select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
system("mkdir -p Analysis/ATAC_seq/Union/R_output_Common_2/NK")
write.table(NK,file="Analysis/ATAC_seq/Union/R_output_Common_2/NK/NK.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#-----FC2 Common to CD4 and CD8, but unique to all others, including NK cells---------------
CD4.CD8 <- reads.X.1 %>% 
  filter((CD4-B)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc) %>% filter((CD4-NK)>fc) %>%
  filter((CD8-B)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% filter((CD8-NK)>fc) %>% select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
system("mkdir -p Analysis/ATAC_seq/Union/R_output_Common_2/CD4_CD8")
write.table(CD4.CD8,file="Analysis/ATAC_seq/Union/R_output_Common_2/CD4_CD8/CD4_CD8.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#-----FC2 Common to NK, CD4, and CD8, but unique to all others---------------
CD4.CD8.NK <- reads.X.1 %>% 
  filter((CD4-B)>fc) %>% filter((CD4-CMP)>fc) %>% filter((CD4-GMP)>fc) %>% filter((CD4-GN)>fc) %>% filter((CD4-Lsk)>fc) %>% filter((CD4-MEP)>fc) %>% filter((CD4-Mono)>fc)  %>% 
  filter((CD8-B)>fc) %>% filter((CD8-CMP)>fc) %>% filter((CD8-GMP)>fc) %>% filter((CD8-GN)>fc) %>% filter((CD8-Lsk)>fc) %>% filter((CD8-MEP)>fc) %>% filter((CD8-Mono)>fc) %>% 
  filter((NK-B)>fc) %>% filter((NK-CMP)>fc) %>% filter((NK-GMP)>fc) %>% filter((NK-GN)>fc) %>% filter((NK-Lsk)>fc) %>% filter((NK-MEP)>fc) %>% filter((NK-Mono)>fc) %>% select(CHR, Start, End, B, CD4, CD8, CMP, GMP, GN, Lsk, MEP, Mono, NK)
system("mkdir -p Analysis/ATAC_seq/Union/R_output_Common_2/CD4_CD8_NK")
write.table(CD4.CD8.NK,file="Analysis/ATAC_seq/Union/R_output_Common_2/CD4_CD8_NK/CD4_CD8_NK.bed",sep="\t", col.names = F, row.names = F, quote=FALSE)

#----Combine modules into heatmap (not reduced)-----
heatmap <- B
heatmap <- bind_rows(heatmap, Lsk)
heatmap <- bind_rows(heatmap, NK)
heatmap <- bind_rows(heatmap, B.Lsk)
heatmap <- bind_rows(heatmap, CD4.CD8)
heatmap <- bind_rows(heatmap, CD4.CD8.NK)

#-----PCA on the modules to reduce dimensions-----
B.pca <- prcomp(B[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- B.pca[[3]]
Lsk.pca <- prcomp(Lsk[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, Lsk.pca[[3]])
NK.pca <- prcomp(NK[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, NK.pca[[3]])
B.Lsk.pca <- prcomp(B.Lsk[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, B.Lsk.pca[[3]])
CD4.CD8.pca <- prcomp(CD4.CD8[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, CD4.CD8.pca[[3]])
CD4.CD8.NK.pca <- prcomp(CD4.CD8.NK[,4:13], center = TRUE, scale. = TRUE) 
heatmap.pca <- rbind(heatmap.pca, CD4.CD8.NK.pca[[3]])

#----Plot Heatmap with Reduced Dimensions----
#---Create Labels for Rows----
labels = nrow(get(positiveResults[1]))
for(i in 2:nrow(positiveResults)){
  labels = c(labels, nrow(get(positiveResults[i])))
}
labels
labels.row = labels

pdf('Analysis/ATAC_seq/Union/R_output_Common_2/Heatmap_Test_1_Reduced.pdf', onefile=FALSE)
pheatmap(heatmap.pca, cellwidth = 30, cellheight = 30, scale = "none",
         show_rownames = T, show_colnames = T, 
         main = "ATAC-seq Peaks",
         cluster_rows = FALSE, cluster_cols = FALSE, labels_row = labels.row)
dev.off()
