#wd <- commandArgs(TRUE)
#setwd(wd)
pcr <- read.delim("/mnt/data1/John/Pioneer-Factors/info/pcr_duplicates.txt", header=FALSE, sep=" ")
colnames(pcr) <- c("Seq", "Cell", "PCR_Remove", "Raw")
head(pcr)
prop <- pcr$PCR_Remove/pcr$Raw
pcr <- cbind(pcr, prop)
mean(pcr$prop)
