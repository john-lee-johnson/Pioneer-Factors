##Set working directory
wd <- commandArgs(TRUE)
setwd(wd)
setwd("/mnt/data1/John/Pioneer_Factors")
library(Hmisc)

#---------Normalized to Total Number of Reads (60981718)-----------------------
MNase.DP.Tcf1 <- read.delim("Analysis/MNase_seq/DP/Tcf1_histfile_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Tcf1 <- read.delim("Analysis/MNase_seq/DN/Tcf1_histfile_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.Ets1 <- read.delim("Analysis/MNase_seq/DP/Ets1_histfile_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Ets1 <- read.delim("Analysis/MNase_seq/DN/Ets1_histfile_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.Runx1 <- read.delim("Analysis/MNase_seq/DP/Runx1_histfile_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Runx1 <- read.delim("Analysis/MNase_seq/DN/Runx1_histfile_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.Gata3 <- read.delim("Analysis/MNase_seq/DP/Gata3_histfile_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Gata3 <- read.delim("Analysis/MNase_seq/DN/Gata3_histfile_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.EML <- read.delim("Analysis/MNase_seq/DP/EML_histfile_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.EML <- read.delim("Analysis/MNase_seq/DN/EML_histfile_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

dev.off()
y1=2
y2=10
mypar(2,2)
plot(MNase.DN.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)


plot(MNase.DN.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Ets1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Ets1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

plot(MNase.DN.Runx1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Runx1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Runx1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Runx1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

plot(MNase.DN.Gata3[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Gata3", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Gata3[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Gata3", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

#---------Normalized To Total Number of Reads (60981718) Shifted bam-----------------
MNase.DP.Tcf1 <- read.delim("Analysis/MNase_seq/DP/Tcf1_histfile_shift_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Tcf1 <- read.delim("Analysis/MNase_seq/DN/Tcf1_histfile_shift_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.Ets1 <- read.delim("Analysis/MNase_seq/DP/Ets1_histfile_shift_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Ets1 <- read.delim("Analysis/MNase_seq/DN/Ets1_histfile_shift_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.Runx1 <- read.delim("Analysis/MNase_seq/DP/Runx1_histfile_shift_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Runx1 <- read.delim("Analysis/MNase_seq/DN/Runx1_histfile_shift_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.Gata3 <- read.delim("Analysis/MNase_seq/DP/Gata3_histfile_shift_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Gata3 <- read.delim("Analysis/MNase_seq/DN/Gata3_histfile_shift_reads.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.EML <- read.delim("Analysis/MNase_seq/DP/EML_histfile_shift_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.EML <- read.delim("Analysis/MNase_seq/DN/EML_histfile_shift_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

head(MNase.DP.Tcf1)

dev.off()
y1=2
y2=10
mypar(2,2)
plot(MNase.DN.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)


plot(MNase.DN.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Ets1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Ets1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

plot(MNase.DN.Runx1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Runx1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Runx1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Runx1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

plot(MNase.DN.Gata3[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Gata3", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Gata3[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Gata3", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

#plot(MNase.DN.EML[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Tcf1 EML", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
#par(new = TRUE)
#plot(MNase.DP.EML[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Tcf1 EML", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

#---------Normalized To Local Area raw bam-----------------
MNase.DP.Tcf1 <- read.delim("Analysis/MNase_seq/DP/Tcf1_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Tcf1 <- read.delim("Analysis/MNase_seq/DN/Tcf1_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.Ets1 <- read.delim("Analysis/MNase_seq/DP/Ets1_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Ets1 <- read.delim("Analysis/MNase_seq/DN/Ets1_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.Runx1 <- read.delim("Analysis/MNase_seq/DP/Runx1_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Runx1 <- read.delim("Analysis/MNase_seq/DN/Runx1_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.Gata3 <- read.delim("Analysis/MNase_seq/DP/Gata3_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Gata3 <- read.delim("Analysis/MNase_seq/DN/Gata3_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.EML <- read.delim("Analysis/MNase_seq/DP/EML_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.EML <- read.delim("Analysis/MNase_seq/DN/EML_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

head(MNase.DP.Tcf1)

dev.off()
y1=0.0
y2=0.1
mypar(2,2)
plot(MNase.DN.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)


plot(MNase.DN.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Ets1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Ets1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

plot(MNase.DN.Runx1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Runx1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Runx1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Runx1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

plot(MNase.DN.Gata3[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Gata3", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Gata3[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Gata3", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

#plot(MNase.DN.EML[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Tcf1 EML", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
#par(new = TRUE)
#plot(MNase.DP.EML[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Tcf1 EML", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

#---------Normalized To Local Area shift bam-----------------
MNase.DP.Tcf1 <- read.delim("Analysis/MNase_seq/DP/Tcf1_histfile_shift_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Tcf1 <- read.delim("Analysis/MNase_seq/DN/Tcf1_histfile_shift_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.Ets1 <- read.delim("Analysis/MNase_seq/DP/Ets1_histfile_shift_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Ets1 <- read.delim("Analysis/MNase_seq/DN/Ets1_histfile_shift_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.Runx1 <- read.delim("Analysis/MNase_seq/DP/Runx1_histfile_shift_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Runx1 <- read.delim("Analysis/MNase_seq/DN/Runx1_histfile_shift_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.Gata3 <- read.delim("Analysis/MNase_seq/DP/Gata3_histfile_shift_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Gata3 <- read.delim("Analysis/MNase_seq/DN/Gata3_histfile_shift_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

MNase.DP.EML <- read.delim("Analysis/MNase_seq/DP/EML_histfile_shift_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.EML <- read.delim("Analysis/MNase_seq/DN/EML_histfile_shift_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))


dev.off()
y1=0.0
y2=0.1
mypar(2,2)
plot(MNase.DN.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)


plot(MNase.DN.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Ets1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Ets1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

plot(MNase.DN.Runx1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Runx1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Runx1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Runx1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)

plot(MNase.DN.Gata3[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Gata3", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Gata3[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Gata3", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)





#--------------1E7---------------------------------------
MNase.DP.Tcf1 <- read.delim("Analysis/MNase_seq/DP/Tcf1_histfile_norm.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Tcf1 <- read.delim("Analysis/MNase_seq/DN/Tcf1_histfile_norm.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DP.Ets1 <- read.delim("Analysis/MNase_seq/DP/Ets1_histfile_norm.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Ets1 <- read.delim("Analysis/MNase_seq/DN/Ets1_histfile_norm.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
dev.off()
y1=0
y2=1
plot(MNase.DN.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DN.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "orange", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "green", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)

#--------------Hist------------------------------------
MNase.DP.Tcf1 <- read.delim("Analysis/MNase_seq/DP/Tcf1_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Tcf1 <- read.delim("Analysis/MNase_seq/DN/Tcf1_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DP.Ets1 <- read.delim("Analysis/MNase_seq/DP/Ets1_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Ets1 <- read.delim("Analysis/MNase_seq/DN/Ets1_histfile_hist.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
dev.off()
y1=0.02
y2=0.04
plot(MNase.DN.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DN.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "orange", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "green", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)

#--------------Shift------------------------------------
MNase.DP.Tcf1 <- read.delim("Analysis/MNase_seq/DP/Tcf1_histfile_shift.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Tcf1 <- read.delim("Analysis/MNase_seq/DN/Tcf1_histfile_shift.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DP.Ets1 <- read.delim("Analysis/MNase_seq/DP/Ets1_histfile_shift.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Ets1 <- read.delim("Analysis/MNase_seq/DN/Ets1_histfile_shift.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
dev.off()
y1=5.5
y2=2.5
plot(MNase.DN.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5,xaxt = "n")
par(tcl= -0.2)
axis(1, at = seq(-2000, 2000, by=250), labels=F, lwd=1, lwd.ticks=1)
par(tcl= -0.5)
axis(1, at = seq(-2000, 2000, by=500), labels=seq(-2000, 2000, by=500), lwd=0, lwd.ticks=2)
par(new = TRUE)
plot(MNase.DP.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DN.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "orange", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "green", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)

#--------------CD4 CD8------------------------------------
MNase.DP.Tcf1.CD4.CD8 <- read.delim("Analysis/MNase_seq/DP/CD4_CD8_histfile_shift.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Tcf1.CD4.CD8 <- read.delim("Analysis/MNase_seq/DN/CD4_CD8_histfile_shift.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DP.Tcf1 <- read.delim("Analysis/MNase_seq/DP/Tcf1_histfile_shift.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Tcf1 <- read.delim("Analysis/MNase_seq/DN/Tcf1_histfile_shift.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
dev.off()
y1=5.5
y2=2.5
plot(MNase.DP.Tcf1.CD4.CD8[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "purple", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5,xaxt = "n")
par(tcl= -0.2)
axis(1, at = seq(-2000, 2000, by=250), labels=F, lwd=1, lwd.ticks=1)
par(tcl= -0.5)
axis(1, at = seq(-2000, 2000, by=500), labels=seq(-2000, 2000, by=500), lwd=0, lwd.ticks=2)
par(new = TRUE)
plot(MNase.DN.Tcf1.CD4.CD8[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "black", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DN.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)

plot(MNase.DN.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "orange", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "green", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)

MNase.DP.Runx1 <- read.delim("Analysis/MNase_seq/DP/Runx1_histfile.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
MNase.DN.Runx1 <- read.delim("Analysis/MNase_seq/DN/Runx1_histfile.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))


#MNase.DP.normalize <- read.delim("Analysis/MNase_seq/DP/normalize.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))
#MNase.DN.normalize <- read.delim("Analysis/MNase_seq/DN/normalize.txt", skip=1, sep = "\t", col.names=c("Distance", "Coverage", "Plus", "Minus"))

#MNase.DP[,2] <- MNase.DP[,2]/MNase.DP.normalize[,3]


dev.off()
y1=0.3
y2=1
plot(MNase.DN.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "blue", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Tcf1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "red", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)


#dev.off()
#y1=0.2
#y2=0.8
plot(MNase.DN.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "orange", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
plot(MNase.DP.Ets1[,c("Distance","Coverage")], type="l", xlim=range(-2000,2000), ylim=range(c(y1,y2)), col = "green", main="Nucleosome Occupancy for Tcf1", xlab="Distance From Center", ylab="Nucleosome Coverage",lwd=2.5)
par(new = TRUE)
