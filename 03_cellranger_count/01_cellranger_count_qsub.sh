#!/bin/bash

# path to the cell ranger script
cellranger_script=/gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/scripts/GB-MN-1318/03_cellranger_count/cellranger_count.sh

#make the results output directory
mkdir /gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/results/cellranger_count

# qsub cell ranger count script for all 16 independent libraries
for d in /gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/data/fastq/YH_MN02/*_R1*; do
	smp_filename=$(basename $d)
	smp=$(awk -F'_S[1-9]' '{print $1}' <<< $smp_filename)
	smp_number=$(grep -Po "_S\d{1,2}_" <<<"$smp_filename")
	smp_number="${smp_number#?}" #removes first character
	smp_id="${smp_number}${smp}"
	qsub $cellranger_script $smp_id $smp
done
