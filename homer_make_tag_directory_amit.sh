#!/bin/bash
cd $dir0/Analysis/Homer/Tag_Directories/Amit
batchMakeTagDirectory.pl $paralleldir/homer_key_file.txt -cpu 40 -genome mm10 -format sam