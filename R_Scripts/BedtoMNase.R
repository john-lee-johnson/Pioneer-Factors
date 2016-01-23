setwd("/mnt/data1/John/Pioneer_Factors")

MNase <- read.delim("Data/MNase_seq/DP/DP_MNase_Seq_mm10_GSM1359852.bed", skip=0, sep = "\t", col.names=c("CHR", "Start", "End", "Weird", "Score", "Strand"))
head(MNase)