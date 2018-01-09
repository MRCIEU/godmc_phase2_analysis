#!/bin/bash

#250MB files
tail -n +2 HiC.GM12878.correlations.txt | split -C 250m --numeric-suffixes -a 3 - ~/GoDMC_Analysis/Hi-C/Grubert2015/chunks/grubert_
for file in ~/GoDMC_Analysis/Hi-C/Grubert2015/chunks/grubert_*
do
    head -n 1 HiC.GM12878.correlations.txt > tmp_file
    cat $file >> tmp_file
    mv -f tmp_file $file
	gzip $file
done
