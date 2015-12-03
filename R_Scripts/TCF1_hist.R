#-------------------------------------------------------------------------------------
# Our objective is to find different clusters of TCF binding and do motif analysis based on
# distinct chromatin states and perform motif analysis.
#-------------------------------------------------------------------------------------
rm(list=ls())
library(GenomicRanges)
library(gplots)
library(pheatmap)
library(reshape)
library(plyr)
library(ggplot2)
library(WGCNA);
library(matrixStats)
library(grid)
library(gridExtra)
#-------------------------------------------------------------------------------------
# Parameters that need to be changed.
beta = 30
Hist_Mark = 'Tcf1'
rep_file = 'Thy_Tcf1_macs_peaks.bedgraph'
#-------------------------------------------------------------------------------------
# Required for multi-threading
#enableWGCNAThreads()
#-------------------------------------------------------------------------------------

Read_HistoneData <- function(SE_file)
{
  tag_cout = unique(read.table(paste(filedir,SE_file,sep=""),
                           header=F,stringsAsFactors=F,sep='\t'))
  return(tag_cout[,5])
}
Plot_motif_count <- function(cluster_number)
{
  #input_file=paste('/Cluster',cluster_number,'_TCF1_motif.count',sep='')
  #output_file=paste('MotifHistogram_HomerAnnotatePeaks_Cluster',cluster_number,'.pdf',sep='')
  input_file="/Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8/CD4_CD8_intersect_motif.count"
  output_file="/Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8/CD4_CD8_intersect_motif.pdf"
  Motif_counts = read.table(paste(mydir,input_file,sep=''),
                            header=T,stringsAsFactors=F,sep='\t')
  
  coordinates =Motif_counts[,1]
  mymotif='ETS'
  Motifs1 = data.frame(coordinates, Motif_counts[,2])
  colnames(Motifs1) = c('coord','tpm')
  Motifs1$id = mymotif
  mymotif='RUNX'
  Motifs2 = data.frame(coordinates, Motif_counts[,5])
  colnames(Motifs2) = c('coord','tpm')
  Motifs2$id = mymotif
  mymotif='TCF'
  Motifs3 = data.frame(coordinates,Motif_counts[,8])
  colnames(Motifs3) = c('coord','tpm')
  Motifs3$id = mymotif
  
  mymotif='SP1'
  Motifs4 = data.frame(coordinates,Motif_counts[,11])
  colnames(Motifs4) = c('coord','tpm')
  Motifs4$id = mymotif
  
  mymotif='KLF'
  Motifs5 = data.frame(coordinates,Motif_counts[,14])
  colnames(Motifs5) = c('coord','tpm')
  Motifs5$id = mymotif
  
  myplot =rbind(Motifs1,
                Motifs2,
                Motifs3,
          #      Motifs4,
                Motifs5)
  
  p = ggplot(myplot) +
  geom_line(aes(x=coord, y=tpm, colour=id), size=0.7)  +
  labs(x="distance to TCF binding (bps)", y='Motif Counts') + 
  ggtitle(paste("Motif counts cluster ", cluster_number,sep='')) 

  ggsave(output_file)
  return(p)
}
#-------------------------------------------------------------------------------------
filedir="/mnt/data1/Golnaz/Pioneer_TF/Idea1_Kmeans_Homer/"
wd <- commandArgs(TRUE)
setwd(wd)
#setwd("/mnt/data1/John/Pioneer_Factors")
mydir=getwd()
#-------------------------------------------------------------------------------------
# just to have the coordinates
p1=Plot_motif_count(1)
p2=Plot_motif_count(2)
p3=Plot_motif_count(3)
p4=Plot_motif_count(4)
p5=Plot_motif_count(5)
p6=Plot_motif_count(6)
p7=Plot_motif_count(7)
p8=Plot_motif_count(8)
pdf('MotifCounts_At_TCFBindingClusters.pdf')
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol=2, nrow =4)
dev.off()
