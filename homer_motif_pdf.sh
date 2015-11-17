#!/bin/bash
cd /mnt/data1/John/Projects/Pioneer_Factors/Analysis/ChIP_seq/motif
for i in $(ls */*Results/*.motif); do
cd $(echo $i | rev | cut -d'/' -f2,3,4,5,6,7,8,9,10 | rev)
motif2Logo.pl $(echo $i | rev | cut -d'/' -f1 | rev) -pdf
echo $(echo $i | rev | cut -d'/' -f1 | rev) -pdf
cd /mnt/data1/John/Projects/Pioneer_Factors/Analysis/ChIP_seq/motif
#motif2Logo.pl $i -pdf
done

