#!/bin/bash
cd /mnt/data1/John/pioneer/analysis/chip_seq/tcf1/macs/motif
#for i in `ls -d *10_*`; do
#for i in `ls -d *4_*`; do
for i in `ls -d *2_*`; do
cd $i
echo $i
changeNewLine.pl ${i}.bed
bedtools getfasta -fi /mnt/data0/John/GRCm38p4_mm10_pa_only.fa  -bed ${i}.bed -fo peaks.fa
meme-chip -oc . -index-name meme-chip.html -meme-p 40 -order 1 -db /mnt/data1/bin/meme_4.10.2/db/motif_databases/MOUSE/uniprobe_mouse.meme -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 10 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-local -centrimo-ethresh 10.0 peaks.fa
cd ..
done
