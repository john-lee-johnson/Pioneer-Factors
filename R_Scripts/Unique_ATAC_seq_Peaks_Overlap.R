setwd("/mnt/data1/John/Pioneer_Factors")
library(rafalib)
txtpath="Analysis/ATAC_seq/Union/R_output_Common_1.5/"
tf.overlap <- read.delim(paste0(txtpath,"tcf1_ets1.txt"), sep = " ", header=FALSE)

pdf('Binding_Proportions_TF_1.5FC.pdf', width=18, height=11)
mypar(3,5)
for (i in 1:nrow(tf.overlap)){
  assign(paste0("t",i),tf.overlap[i,2])
  total <- eval(parse(text = paste0("t",i)))
  assign(paste0("p",i),c(tf.overlap[i,4]/total,tf.overlap[i,6]/total,tf.overlap[i,8]/total)) 
  assign(paste0("m",i),tf.overlap[i,1]) 
  title = eval(parse(text = paste0("m",i)))
  xlabel = paste0("Total ATAC Peaks= ",eval(parse(text = paste0("t",i))))
  barplot(eval(parse(text = paste0("p",i))), main=title, xlab=xlabel, ylab="Proportion binding at ATAC-seq sites", names.arg=c("ETS1_DN", "ETS1_DP", "TCF1"), ylim=c(0,0.8))
}  
dev.off()

txtpath="Analysis/ATAC_seq/Union/R_output_Common_Merge1.5/"
tf.overlap <- read.delim(paste0(txtpath,"tcf1_ets1.txt"), sep = " ", header=FALSE)
nrow(tf.overlap)
pdf('Binding_Proportions_TF_1.5FC_Merge.pdf', width=18, height=11)
mypar(3,5)
for (i in 1:nrow(tf.overlap)){
  assign(paste0("t",i),tf.overlap[i,2])
  total <- eval(parse(text = paste0("t",i)))
  assign(paste0("p",i),c(tf.overlap[i,4]/total,tf.overlap[i,6]/total,tf.overlap[i,8]/total)) 
  assign(paste0("m",i),tf.overlap[i,1]) 
  title = eval(parse(text = paste0("m",i)))
  xlabel = paste0("Total ATAC Peaks= ",eval(parse(text = paste0("t",i))))
  barplot(eval(parse(text = paste0("p",i))), main=title, xlab=xlabel, ylab="Proportion binding at ATAC-seq sites", names.arg=c("ETS1_DN", "ETS1_DP", "TCF1"), ylim=c(0,0.8))
}  
dev.off()

txtpath="Analysis/ATAC_seq/Overlap/"
tf.overlap.cluster <- read.delim(paste0(txtpath,"atac_chip_overlap_clusters_0.txt"), sep = " ", header=FALSE, col.names=c("Cluster","totalCluster","Cell", "ChIP", "Investigator", "Bound", "Total","Proportion"))
head(tf.overlap.cluster)
dp <- tf.overlap.cluster[tf.overlap.cluster$Cell=="DP" | tf.overlap.cluster$Cell=="Thy",]
dp <- dp[dp$ChIP!="PU.1",] 
p1 <- ggplot(dp, aes(x=Cluster, y=Proportion)) + geom_point(aes(color=ChIP),shape = 1, size = 4) + theme_bw() + theme(text = element_text(size=20))
p1
p2 <- p1 + scale_x_discrete(labels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
p2
p3 <- p2 +
  geom_point(aes(color=ChIP),size = 4.5, shape = 1) +
  geom_point(aes(color=ChIP),size = 4, shape = 1) +
  geom_point(aes(color=ChIP),size = 3.5, shape = 1)
geom_point(aes(color=ChIP),size = 3, shape = 1)
pdf(paste0(txtpath,'Binding_Proportions_TF_0_Cluster15.pdf'), width=11, height=8.5)
postscript(paste0(txtpath,'Binding_Proportions_TF_0_Cluster15.eps'), width=11, height=8.5)
p3
dev.off()

txtpath="Analysis/ATAC_seq/Overlap/"
tf.overlap.cluster <- read.delim(paste0(txtpath,"atac_chip_overlap_clusters_0_merge.txt"), sep = " ", header=FALSE, col.names=c("Cluster","totalCluster","Cell", "ChIP", "Investigator", "Bound", "Total","Proportion","Average"))
head(tf.overlap.cluster)
tcf1 <- tf.overlap.cluster %>% filter(Cell=="Thy") %>% filter(ChIP=="Tcf1")
tcf1$Zscore <- scale(tcf1$Average, center = TRUE, scale = TRUE)
tcf1
tf.overlap.cluster[tf.overlap.cluster$ChIP=="Tcf1",]
dp <- tf.overlap.cluster[tf.overlap.cluster$Cell=="DP" | tf.overlap.cluster$Cell=="Thy" |  tf.overlap.cluster$Cell=="EML",]
dp <- dp[dp$ChIP!="PU.1",] 
p1 <- ggplot(dp, aes(x=Cluster, y=Proportion)) + geom_point(aes(color=Average,shape = ChIP), size = 3) + theme_bw() + theme(text = element_text(size=20))
p1
p2 <- p1 + scale_x_discrete(labels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
p2
p3 <- p2 +
  geom_point(aes(color=ChIP),size = 4.5, shape = 1) +
  geom_point(aes(color=ChIP),size = 4, shape = 1) +
  geom_point(aes(color=ChIP),size = 3.5, shape = 1)
geom_point(aes(color=ChIP),size = 3, shape = 1)
postscript(paste0(txtpath,'Binding_Proportions_TF_0_Cluster15_merge.eps'), width=11, height=8.5)
p3
dev.off()

dp

clusterTotal <- tf.overlap.cluster[tf.overlap.cluster$ChIP=="Tcf1",]
clusterTotal <- clusterTotal[clusterTotal$Cell!="EML",]
sum(clusterTotal$Bound)
