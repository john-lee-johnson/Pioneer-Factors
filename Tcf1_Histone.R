library(NMF)
library(WGCNA)
library(RColorBrewer)
library(rtracklayer)
library(biomaRt)
library(Rsamtools)
library(refGenome)
library(dynamicTreeCut)
library(ChIPseeker)
library(ChIPpeakAnno)
library(pheatmap)
#setwd("~/Desktop/John/pioneer")
#setwd("~/Desktop/John/pioneer")
setwd("/mnt/data1/John/pioneer")
#setwd("~/Desktop/John/Pioneer-Factors/files")
#setwd("/mnt/data1/John/Pioneer-Factors/files")

mapped.reads=read.delim("../Pioneer-Factors/files/mapped_reads_histone.txt", sep=" ",header=FALSE)
colnames(mapped.reads) <- c("Cell", "Mark", "Reads")
bedPath = "analysis/chip_seq/tcf1/tag/"
myFiles <- list.files(path=bedPath, pattern="*tag.bed")
cell <- NULL
chr <- as.character(read.table(paste0(bedPath,myFiles[1]),header=FALSE)[,1])
for (i in myFiles ) {
  cell <- unlist(strsplit(i,"_")[1])[1]
  mark <- unlist(strsplit(i,"_")[1])[2]
  name <- paste0(cell, "_", mark)
  file = paste0(bedPath, name, "_Tcf1_Thy_tag.bed")
  norm <- mapped.reads$Reads[mapped.reads$Cell==cell & mapped.reads$Mark==mark]
  data <- read.delim(file, sep="\t", header=FALSE)
  colnames(data) <- c("CHR", "Start", "End", "score")
  data$score <- data$score/norm
  write.table(data, file=paste0(bedPath, name, "_", "Tcf1_Thy_tag_norm.bed"), sep="\t", quote= FALSE, row.names=FALSE, col.names=FALSE)
  data = NULL
  #data <- import.bedGraph(file, trackLine=FALSE, genome="mm10")
  data <- read.delim(file = paste0(bedPath, name, "_", "Tcf1_Thy_tag_norm.bed"),header=FALSE)
  #data$score <- data$score/norm
  colnames(data) <- c("CHR","Start","End","score")
  data$score <- log2(data$score+1)
  write.table(data, file=paste0(bedPath, name, "_", "Tcf1_Thy_tag_norm_log.bed"), sep="\t", quote= FALSE, row.names=FALSE, col.names=FALSE)
  assign(name,data)
}
Tcf1 <- read.delim("analysis/chip_seq/tcf1/macs/Thy_Tcf1_macs_peaks.bed", sep="\t", header=FALSE)
colnames(Tcf1) <- c("CHR", "Start", "End", "Peak", "Score")
Tcf1.sort <- Tcf1[order(Tcf1$Score, decreasing = TRUE),]
head(Tcf1.sort)
Tcf1.sort["Strand"] <- 0
num=trunc(nrow(Tcf1.sort)/10)
Tcf1.sort.window <- Tcf1.sort
Tcf1.sort.window$Start <- trunc((Tcf1.sort$Start + Tcf1.sort$End)/2-100)
Tcf1.sort.window$End <- trunc((Tcf1.sort$Start + Tcf1.sort$End)/2+100)
start=1; end=num
for (i in 1:10) {
  system(paste0("mkdir -p analysis/chip_seq/tcf1/macs/motif/Thy_Tcf1_macs_peaks_10_", i))
  write.table(Tcf1.sort.window[start:end,], file=paste0("analysis/chip_seq/tcf1/macs/motif/Thy_Tcf1_macs_peaks_10_", i, "/Thy_Tcf1_macs_peaks_10_", i , ".bed"), sep="\t", quote= FALSE, row.names=FALSE, col.names=FALSE)
  start=start+num
  end=end+num
}
num=trunc(nrow(Tcf1.sort)/4)
Tcf1.sort.window <- Tcf1.sort
Tcf1.sort.window$Start <- trunc((Tcf1.sort$Start + Tcf1.sort$End)/2-100)
Tcf1.sort.window$End <- trunc((Tcf1.sort$Start + Tcf1.sort$End)/2+100)
start=1; end=num
for (i in 1:4) {
  system(paste0("mkdir -p analysis/chip_seq/tcf1/macs/motif/Thy_Tcf1_macs_peaks_4_", i))
  write.table(Tcf1.sort.window[start:end,], file=paste0("analysis/chip_seq/tcf1/macs/motif/Thy_Tcf1_macs_peaks_4_", i, "/Thy_Tcf1_macs_peaks_4_", i , ".bed"), sep="\t", quote= FALSE, row.names=FALSE, col.names=FALSE)
  start=start+num
  end=end+num
}
num=trunc(nrow(Tcf1.sort)/2)
Tcf1.sort.window <- Tcf1.sort
Tcf1.sort.window$Start <- trunc((Tcf1.sort$Start + Tcf1.sort$End)/2-100)
Tcf1.sort.window$End <- trunc((Tcf1.sort$Start + Tcf1.sort$End)/2+100)
start=1; end=num
for (i in 1:2) {
  system(paste0("mkdir -p analysis/chip_seq/tcf1/macs/motif/Thy_Tcf1_macs_peaks_2_", i))
  write.table(Tcf1.sort.window[start:end,], file=paste0("analysis/chip_seq/tcf1/macs/motif/Thy_Tcf1_macs_peaks_2_", i, "/Thy_Tcf1_macs_peaks_2_", i , ".bed"), sep="\t", quote= FALSE, row.names=FALSE, col.names=FALSE)
  start=start+num
  end=end+num
}


#Depreciated Gene Annotation
#------------------------------------------------------------------------------------------------------------------------
#Tcf1.bedgraph <- read.delim("analysis/chip_seq/tcf1/macs/Thy_Tcf1_macs_peaks.bedgraph", sep="\t", header=FALSE)
#colnames(Tcf1.bedgraph) <- c("CHR", "Start", "End", "Score")
#Tcf1.sort.bedgraph <- Tcf1.bedgraph[order(Tcf1.bedgraph$Score, decreasing = TRUE),]
#Tcf1_annotate <- read.delim("analysis/chip_seq/tcf1/tag/Thy_Tcf1_window_sort_homer_annotate_sort.bed", sep="\t", header=TRUE)
#Tcf1_annotate_genes <- Tcf1_annotate[,c("Chr", "Start", "End", "Gene.Name", "Gene.Description")]
#Tcf1_genes <- GRanges(seqnames = Rle(Tcf1_annotate_genes$Chr), ranges = IRanges(start = Tcf1_annotate_genes$Start, end = Tcf1_annotate_genes$End), names = Tcf1_annotate_genes$Gene.Name, description = Tcf1_annotate_genes$Gene.Description)
#ensembl <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
#genes<-getBM(attributes=c("external_gene_name","chromosome_name","start_position","end_position"),mart=ensembl)
#genes.clean <- genes[ which( genes$chromosome_name=="1" | genes$chromosome_name=="2" | genes$chromosome_name=="3" | genes$chromosome_name=="4" | genes$chromosome_name=="5" | genes$chromosome_name=="6" | genes$chromosome_name=="7" | genes$chromosome_name=="8" | genes$chromosome_name=="9" | genes$chromosome_name=="10" | genes$chromosome_name=="11" | genes$chromosome_name=="12" | genes$chromosome_name=="13" | genes$chromosome_name=="14" | genes$chromosome_name=="15" | genes$chromosome_name=="16" | genes$chromosome_name=="17" | genes$chromosome_name=="18" | genes$chromosome_name=="19" | genes$chromosome_name=="X" | genes$chromosome_name=="Y" | genes$chromosome_name=="MT") , ]; genes.clean$chromosome_name[which(genes.clean$chromosome_name=="MT")]<-"M"; genes.clean<- genes.clean[order(genes.clean$chromosome_name, decreasing=FALSE),]
#genes_int <- GRanges(seqnames = Rle(paste0("chr",genes.clean$chromosome_name)), ranges = IRanges(start = genes.clean$start_position, end = genes.clean$end_position), names = genes.clean$external_gene_name, genome="mm10")
#overlaps <- findOverlaps(LTHSC_H3K4me1, genes_int)
#LTHSC_H3K4me1$names<-NA
#for (i in 1:length(overlaps) ) {
#  a <- queryHits(overlaps)[i]
#  b <- subjectHits(overlaps)[i]
#  LTHSC_H3K4me1$names[a]<- genes_int$names[b]
#}
#overlaps_Tcf1 <- findOverlaps(LTHSC_H3K4me1, Tcf1_genes)
#LTHSC_H3K4me1$names<-NA
#for (i in 1:length(overlaps_Tcf1) ) {
#  a <- queryHits(overlaps_Tcf1)[i]
#  b <- subjectHits(overlaps_Tcf1)[i]
#  LTHSC_H3K4me1$names[a]<- Tcf1_genes$names[b]
#  LTHSC_H3K4me1$description[a]<- Tcf1_genes$description[b]
#}
#rownames(dat.me1)<-LTHSC_H3K4me1$names

hmcol <- colorRampPalette(brewer.pal(6, "RdBu"))(10)
hmcol <- rev(hmcol)

dat.me1 <- cbind(LTHSC_H3K4me1$score, STHSC_H3K4me1$score, MPP_H3K4me1$score, CLP_H3K4me1$score, CD4_H3K4me1$score, CD8_H3K4me1$score); colnames(dat.me1)<-c("LTHSC", "STHSC", "MPP", "CLP", "CD4", "CD8"); ind = NULL; ind <- apply(dat.me1, 1, var) == 0; dat.me1 <- dat.me1[!ind,]
dat.me1.bedgraph <- cbind(as.character(LTHSC_H3K4me1$CHR), LTHSC_H3K4me1$Start, LTHSC_H3K4me1$End, LTHSC_H3K4me1$score, STHSC_H3K4me1$score, MPP_H3K4me1$score, CLP_H3K4me1$score, CD4_H3K4me1$score, CD8_H3K4me1$score); colnames(dat.me1.bedgraph)<-c("CHR", "Start", "End", "LTHSC", "STHSC", "MPP", "CLP", "CD4", "CD8"); dat.me1.bedgraph <- dat.me1.bedgraph[!ind,]

dissim.me1 = as.dist(1 - cor(t(dat.me1)))
dendro.me1 <- hclust(dissim.me1, method="average")


DetectedColors.me1 = NULL;

DetectedColors.me1 = cbind(DetectedColors.me1, labels2colors(cutreeDynamic(dendro.me1, cutHeight = 1.5, minClusterSize = 200, method = "hybrid", deepSplit = FALSE, pamStage = TRUE,  distM = as.matrix(dissim.me1), maxPamDist = 0, verbose = 0)));

Methods.me1 = c("Dynamic Hybrid");

plot(dendro.me1, main = "Clustering dendrogram and module colors", ylab = "Difference");
plotDendroAndColors(dendro.me1, DetectedColors.me1, groupLabels = Methods.me1, main="");
dev.off()
row.order <- rev(dendro.me1$order)
myColors.me1 <- list(levels(factor(DetectedColors.me1[rev(dendro.me1$order),])))
# following code limits the lowest and highest color to 5%, and 95% of your range, respectively
quantile.range <- quantile(dat.me1, probs = seq(0, 1, 0.01))
quantile.range
palette.breaks <- seq(quantile.range["5%"], quantile.range["99%"], 0.1)
palette.breaks

#breaks for the core of the distribution
breaks=seq(0, 4, by=0.2) #41 values
#now add outliers
color.palette  <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(length(breaks)-1))

pheatmap(dat.me1[row.order,], color = color.palette, breaks=breaks, cellwidth=50, scale = "none", cluster_rows=FALSE, cluster_cols=FALSE, legend=TRUE, main=paste("H3K4me1 marks + nRow = ", nrow(dat.me1)),filename="/mnt/data1/John/heatmap_me1_no_scale.png", labels_row="")
pheatmap(dat.me1[row.order,], color = color.palette, breaks=breaks, cellwidth=50, scale = "none", cluster_rows=TRUE, cluster_cols=FALSE, legend=TRUE, main=paste("H3K4me1 marks + nRow = ", nrow(dat.me1)),filename="/mnt/data1/John/heatmap_me1_no_scale_cluster.png", labels_row="")
breaks=seq(0, 6, by=0.2) #41 values
color.palette  <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(length(breaks)-1))
pheatmap(dat.me2[row.order,], color = color.palette, breaks=breaks, cellwidth=50, scale = "none", cluster_rows=FALSE, cluster_cols=FALSE, legend=TRUE, main=paste("H3K4me2 marks + nRow = ", nrow(dat.me2)),filename="/mnt/data1/John/heatmap_me2_no_scale.png", labels_row="")
pheatmap(dat.me2[row.order,], color = color.palette, breaks=breaks, cellwidth=50, scale = "none", cluster_rows=TRUE, cluster_cols=FALSE, legend=TRUE, main=paste("H3K4me2 marks + nRow = ", nrow(dat.me2)),filename="/mnt/data1/John/heatmap_me2_no_scale_cluster.png", labels_row="")

aheatmap(dat.me1[row.order,], color = hmcol, cellwidth=50, scale = "row", annRow=DetectedColors.me1[row.order,], annColors=myColors.me1, Rowv=NA, Colv=NA, legend=TRUE, main=paste("H3K4me1 marks + nRow = ", nrow(dat.me1)),filename="/mnt/data1/John/heatmap_me1_scale.png")
aheatmap(dat.me2[row.order,], color = hmcol, cellwidth=50, scale = "row", annRow=DetectedColors.me1[row.order,], annColors=myColors.me1, Rowv=NA, Colv=NA, legend=TRUE, main=paste("H3K4me2 marks + nRow = ", nrow(dat.me2)),filename="/mnt/data1/John/heatmap_me2_scale.png")
write.table(dat.me1.bedgraph[row.order,],file="analysis/chip_seq/tcf1/R_output/H3K4me1_Heatmap.bedgraph", sep="\t", col.names = T, row.names = F, quote=FALSE)
write.table(dat.me2.bedgraph[row.order,],file="analysis/chip_seq/tcf1/R_output/H3K4me2_Heatmap.bedgraph", sep="\t", col.names = T, row.names = F, quote=FALSE)
head(dat.me2[row.order,])
yellow.genes <- unique(rownames(dat.me1[which(DetectedColors=="yellow"),]))
yellow.genes <- yellow.genes[!is.na(yellow.genes)]
write.table(yellow.genes,"/mnt/data1/John/yellowgenes.txt",col.names = F, row.names = F, quote=FALSE)

blue.genes <- unique(rownames(dat.me1[which(DetectedColors=="blue"),]))
blue.genes <- blue.genes[!is.na(blue.genes)]
write.table(blue.genes,"/mnt/data1/John/bluegenes.txt",col.names = F, row.names = F, quote=FALSE)

brown.genes <- unique(rownames(dat.me1[which(DetectedColors=="brown"),]))
brown.genes <- brown.genes[!is.na(brown.genes)]
write.table(brown.genes,"/mnt/data1/John/browngenes.txt",col.names = F, row.names = F, quote=FALSE)

rownames(dat.me1)
turquoise.genes <- unique(rownames(dat.me1[which(DetectedColors.me1=="turquoise"),]))
turquoise.genes <- turquoise.genes[!is.na(turquoise.genes)]
write.table(turquoise.genes,"/mnt/data1/John/turquoisegenes.txt",col.names = F, row.names = F, quote=FALSE)

#--------------H3K4me2----------------------

dat.me2 <- cbind(LTHSC_H3K4me2$score, STHSC_H3K4me2$score, MPP_H3K4me2$score, CLP_H3K4me2$score, CD4_H3K4me2$score, CD8_H3K4me2$score); colnames(dat.me2)<-c("LTHSC", "STHSC", "MPP", "CLP", "CD4", "CD8"); dat.me2 <- dat.me2[!ind,]
dat.me2.bedgraph <- cbind(as.character(LTHSC_H3K4me2$CHR), LTHSC_H3K4me2$Start, LTHSC_H3K4me2$End, LTHSC_H3K4me2$score, STHSC_H3K4me2$score, MPP_H3K4me2$score, CLP_H3K4me2$score, CD4_H3K4me2$score, CD8_H3K4me2$score); colnames(dat.me2.bedgraph)<-c("CHR", "Start", "End", "LTHSC", "STHSC", "MPP", "CLP", "CD4", "CD8"); dat.me2.bedgraph <- dat.me2.bedgraph[!ind,]

aheatmap(dat.me2[row.order,], color = hmcol, cellwidth=50, scale = "none", annRow=DetectedColors.me1[row.order,], annColors=myColors.me1, Rowv=NA, Colv=NA, legend=TRUE, main=paste("H3K4me2 marks + nRow = ", nrow(dat.me2)),filename="/mnt/data1/John/heatmap_me2.png")
dissim.turq <- 1 - cor(t(turq.me2))
dendro.turq <- hclust(as.dist(dissim.turq))

#--------------H3K4me1 and H3K4me2----------------------------
me1.me2 <- dat.me1+dat.me2
rownames(me1.me2) <- rownames(dat.me1)
dissim.me1.me2 = as.dist(1 - cor(t(me1.me2)))
dendro.me1.me2 <- hclust(dissim.me1.me2, method="average")

DetectedColors.me1.me2 = NULL;

DetectedColors.me1.me2 = cbind(DetectedColors.me1.me2, labels2colors(cutreeDynamic(dendro.me1.me2, cutHeight = 1.5, minClusterSize = 200, method = "hybrid", deepSplit = FALSE, pamStage = TRUE,  distM = as.matrix(dissim.me1.me2), maxPamDist = 0, verbose = 0)));

Methods.me1.me2 = c("Dynamic Hybrid");

plot(dendro.me1.me2, main = "Clustering dendrogram and module colors", ylab = "Difference");
plotDendroAndColors(dendro.me1.me2, DetectedColors.me1.me2, groupLabels = Methods.me1.me2, main="");
row.order.me1.me2 <- rev(dendro.me1.me2$order)
myColors.me1.me2 <- list(levels(factor(DetectedColors.me1.me2[rev(dendro.me1.me2$order),])))
aheatmap(me1.me2[row.order.me1.me2,], color = hmcol, cellwidth=50, scale = "none", annRow=DetectedColors.me1.me2[row.order.me1.me2,], annColors=myColors.me1.me2, Rowv=NA, Colv=NA, legend=TRUE, main=paste("H3K4me1 marks+H3K4me2 marks + nRow = ", nrow(dat.me1)),filename="/mnt/data1/John/heatmap_me1_me2_add.png")

turquoise.genes <- unique(rownames(me1.me2[which(DetectedColors.me1.me2=="turquoise"),]))
turquoise.genes <- turquoise.genes[!is.na(turquoise.genes)]
write.table(turquoise.genes,"/mnt/data1/John/turquoisegenes.txt",col.names = F, row.names = F, quote=FALSE)

brown.genes <- unique(rownames(me1.me2[which(DetectedColors.me1.me2=="brown"),]))
brown.genes <- brown.genes[!is.na(brown.genes)]
write.table(brown.genes,"/mnt/data1/John/browngenes.txt",col.names = F, row.names = F, quote=FALSE)

yellow.genes <- unique(rownames(me1.me2[which(DetectedColors.me1.me2=="yellow"),]))
yellow.genes <- yellow.genes[!is.na(yellow.genes)]
write.table(yellow.genes,"/mnt/data1/John/yellowgenes.txt",col.names = F, row.names = F, quote=FALSE)
#--------------------------H3K4me3-----------------------------
dat.me3 <- cbind(LTHSC_H3K4me3$score, STHSC_H3K4me3$score, MPP_H3K4me3$score, CLP_H3K4me3$score, CD4_H3K4me3$score, CD8_H3K4me3$score); colnames(dat.me3)<-c("LTHSC", "STHSC", "MPP", "CLP", "CD4", "CD8"); dat.me3 <- dat.me3[!ind,]

rownames(dat.me3) <- rownames(dat.me1)
dissim.me3 = as.dist(1 - cor(t(dat.me3)))
dendro.me3 <- hclust(dissim.me3, method="average")

DetectedColors.me3 = NULL;

DetectedColors.me3 = cbind(DetectedColors.me3, labels2colors(cutreeDynamic(dendro.me3, cutHeight = 1.5, minClusterSize = 200, method = "hybrid", deepSplit = FALSE, pamStage = TRUE,  distM = as.matrix(dissim.me3), maxPamDist = 0, verbose = 0)));

Methods.me3 = c("Dynamic Hybrid");

plot(dendro.me3, main = "Clustering dendrogram and module colors", ylab = "Difference");
plotDendroAndColors(dendro.me3, DetectedColors.me3, groupLabels = Methods.me3, main="");
row.order.me3 <- rev(dendro.me3$order)
myColors.me3 <- list(levels(factor(DetectedColors.me3[rev(dendro.me3$order),])))
aheatmap(dat.me3[row.order.me3,], color = hmcol, cellwidth=50, scale = "row", annRow=DetectedColors.me3[row.order.me3,], annColors=myColors.me3, Rowv=NA, Colv=NA, legend=TRUE, main=paste("H3K4me3 marks + nRow = ", nrow(dat.me3)),filename="/mnt/data1/John/heatmap_me3.png")

#--------------------------H3K27Ac-----------------------------
dat.ac <- cbind(LTHSC_H3K27Ac$score, STHSC_H3K27Ac$score, MPP_H3K27Ac$score, CLP_H3K27Ac$score, CD4_H3K27Ac$score, CD8_H3K27Ac$score); colnames(dat.ac)<-c("LTHSC", "STHSC", "MPP", "CLP", "CD4", "CD8"); dat.ac <- dat.ac[!ind,]

rownames(dat.ac) <- rownames(dat.me1)
dissim.ac = as.dist(1 - cor(t(dat.ac)))
dendro.ac <- hclust(dissim.ac, method="average")

DetectedColors.ac = NULL;

DetectedColors.ac = cbind(DetectedColors.ac, labels2colors(cutreeDynamic(dendro.ac, cutHeight = 1.5, minClusterSize = 200, method = "hybrid", deepSplit = FALSE, pamStage = TRUE,  distM = as.matrix(dissim.ac), maxPamDist = 0, verbose = 0)));

Methods.ac = c("Dynamic Hybrid");

plot(dendro.ac, main = "Clustering dendrogram and module colors", ylab = "Difference");
plotDendroAndColors(dendro.ac, DetectedColors.ac, groupLabels = Methods.ac, main="");
row.order.ac <- rev(dendro.ac$order)
myColors.ac <- list(levels(factor(DetectedColors.ac[rev(dendro.ac$order),])))
aheatmap(dat.ac[row.order.ac,], color = hmcol, cellwidth=50, scale = "row", annRow=DetectedColors.ac[row.order.ac,], annColors=myColors.ac, Rowv=NA, Colv=NA, legend=TRUE, main=paste("H3K27Ac marks + nRow = ", nrow(dat.ac)),filename="/mnt/data1/John/heatmap_ac.png")


which(genes_int$names=="Nr4a1") #Nur77
which(LTHSC.me1$names=="Bcl2") #Cell death
which(LTHSC.me1$names=="Bcl11b") #TF for T cell development
which(LTHSC.me1$names=="Gfi1") #Homeostasis of HSC
which(LTHSC.me1$names=="Cd8a")
which(LTHSC.me1$names=="Cd4")
which(LTHSC.me1$names=="Foxp1")
which(LTHSC.me1$names=="Tbx21")
which(LTHSC.me1$names=="Cd40lg")
which(LTHSC.me1$names=="Cd28")
which(LTHSC.me1$names=="Cd3g")
which(LTHSC.me1$names=="Il21r")
which(LTHSC.me1$names=="Cxcr4")
which(LTHSC.me1$names=="Il7r")
which(LTHSC.me1$names=="Nr4a1") #Cell death
which(LTHSC.me1$names=="Gata3") #TF for T cell development
which(LTHSC.me1$names=="Bcl11b")
which(LTHSC.me1$names=="Notch1") #TF for T cell development
which(LTHSC.me1$names=="Tcf7") #TF for T cell development
which(LTHSC.me1$names=="Ccr7")
which(LTHSC.me1$names=="Rag1")
lymphoid.genes <- c("Bcl2","Cd8a","Cd4","Tbx21","Cd40lg","Cd28","Cd3g","Il21r","Cd69","Gata3","Bcl11b","Notch1","Ccr7","Rag1","Cxcr4","Il7r")

