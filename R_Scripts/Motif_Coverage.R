#-------------------------------------------------------------------------------------
# This script will take the motif output from homer and plot the coverage 
#-------------------------------------------------------------------------------------
library(gplots)
library(reshape)
library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)
#-------------------------------------------------------------------------------------
Plot_motif_count <- function(i)
{
  if ( i == 1 ){ 
    input_file="CD4_CD8_intersect_motif.count"
    output_file="CD4_CD8_intersect_motif.pdf"
    motif="CD4 and CD8 Old"
  } else if ( i == 2) { 
    input_file="CD4_CD8_intersect_motif_new.count"
    output_file="CD4_CD8_intersect_motif_new.pdf"
    motif="CD4 and CD8 New"
  } else if ( i == 3) { 
    input_file="NK_CD4_CD8_intersect_motif.count"
    output_file="NK_CD4_CD8_intersect_motif.pdf"
    motif="NK, CD4 and CD8 Old"
  } else if ( i == 4) { 
    input_file="NK_CD4_CD8_intersect_motif_new.count"
    output_file="NK_CD4_CD8_intersect_motif_new.pdf"
    motif="NK, CD4 and CD8 New"
} else if ( i == 5) { 
  input_file="CD4_CD8_intersect_motif_new_chip.count"
  output_file="CD4_CD8_intersect_motif_new_chip.pdf"
  motif="CD4 and CD8 New, Tcf ChIP"
} else if ( i == 6) { 
  input_file="NK_CD4_CD8_intersect_motif_new_chip.count"
  output_file="NK_CD4_CD8_intersect_motif_new_chip.pdf"
  motif="NK, CD4 and CD8 New, Tcf ChIP"
}
  
  Motif_counts = read.table(paste0(mydir,input_file),
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
  
  myplot =rbind(Motifs1,
                Motifs2,
                Motifs3)
  
  p = ggplot(myplot) +
    geom_line(aes(x=coord, y=tpm, colour=id), size=0.5)  +
    labs(x="distance around center (bps)", y='Motif Counts') + 
    ggtitle(paste0(motif," Motifs")) 
  
  ggsave(output_file)
  return(p)
  
}

#-------------------------------------------------------------------------------------
filedir="/mnt/data1/John/Pioneer_Factors/Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8/"
mydir="/mnt/data1/John/Pioneer_Factors/Analysis/ATAC_seq/Union/R_output_Common/CD4_CD8/"
setwd(mydir)
#-------------------------------------------------------------------------------------
# just to have the coordinates
p1=Plot_motif_count(1)
p2=Plot_motif_count(2)
p5=Plot_motif_count(5)

#-------------------------------------------------------------------------------------
filedir="/mnt/data1/John/Pioneer_Factors/Analysis/ATAC_seq/Union/R_output_Common/NK_CD4_CD8/"
mydir="/mnt/data1/John/Pioneer_Factors/Analysis/ATAC_seq/Union/R_output_Common/NK_CD4_CD8/"
setwd(mydir)
#-------------------------------------------------------------------------------------
# just to have the coordinates
p3=Plot_motif_count(3)
p4=Plot_motif_count(4)
p6=Plot_motif_count(6)

#-------------------------------------------------------------------------------------
filedir="/mnt/data1/John/Pioneer_Factors/Analysis/ATAC_seq/Union/R_output_Common/"
mydir="/mnt/data1/John/Pioneer_Factors/Analysis/ATAC_seq/Union/R_output_Common/"
setwd(mydir)

pdf('MotifCounts_At_Common_Peaks.pdf', width=18, height=11)
grid.arrange(p1,p2,p5,p3,p4,p6, ncol=3, nrow =2)

dev.off()

mydir="/mnt/data1/John/Pioneer_Factors/Analysis/ATAC_seq/Union/R_output_Common/Intersect"
setwd(mydir)
# Simple Pie Chart
all <- nrow(read.delim(file="Common_Peaks_All.bed", header=FALSE))
all.Tcf <- nrow(read.delim(file="Common_Peaks_All_Tcf1.bed", header=FALSE))
all.Tcf.window <- nrow(read.delim(file="Common_Peaks_All_Tcf1_Window.bed", header=FALSE))

all
all.Tcf
all.Tcf.window

cd4.cd8 <- nrow(read.delim(file="CD4_CD8_intersect_5_column_merge.bed", header=FALSE))
nk.cd4.cd8 <- nrow(read.delim(file="NK_CD4_CD8_intersect_5_column_merge.bed", header=FALSE))
cd4.cd8.tcf1 <- nrow(read.delim(file="CD4_CD8_Tcf1.bed", header=FALSE))
cd4.cd8
cd4.cd8.tcf1
nk.cd4.cd8.tcf1 <- nrow(read.delim(file="NK_CD4_CD8_Tcf1.bed", header=FALSE))
nk.cd4.cd8.tcf1

cd4.cd8.tcf1.window <- nrow(read.delim(file="CD4_CD8_Tcf1_Window.bed", header=FALSE))
cd4.cd8.tcf1.window
nk.cd4.cd8.tcf1.window <- nrow(read.delim(file="NK_CD4_CD8_Tcf1_Window.bed", header=FALSE))
nk.cd4.cd8.tcf1.window

slices1 <- c(1-all.Tcf/all,all.Tcf/all)
lbls1 <- c("Total Common Peaks","ATAC-seq Peaks with Tcf1")
pie1 <- pie(slices1, labels = lbls1, main="NK, CD4, and CD8 ATAC-seq Peaks that have Tcf1 binding")

slices2 <- c(1-all.Tcf.window/all,all.Tcf.window/all)
lbls2 <- c("Total Common Peaks","ATAC-seq Peaks with Tcf1")
pie2 <- pie(slices2, labels = lbls2, main="NK, CD4, and CD8ATAC-seq Peaks that have Tcf1 binding within 2kb")

pdf('Tcf1_binding_ATAC_seq_peaks.pdf', onefile = FALSE)
mypar(2,1)
pie1 <- pie(slices1, labels = lbls1, main="NK, CD4, and CD8 ATAC-seq Peaks that have Tcf1 binding")
pie2 <- pie(slices2, labels = lbls2, main="NK, CD4, and CD8 ATAC-seq Peaks that have Tcf1 binding within 2kb")
dev.off()

