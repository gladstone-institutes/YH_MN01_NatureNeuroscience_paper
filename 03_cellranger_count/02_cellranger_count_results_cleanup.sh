#!/bin/bash

##########################################################################
## This script takes two command line arguments (in the below order)
## 1. CellRanger results folder path to be cleaned up
## 2. Temporary or trash folder path
##
## Example run:
## 02_cellranger_count_results_cleanup.sh \
## /gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/results/cellranger_count \
## /gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/tmp/cellranger_results
##########################################################################

#check the command line arguments
if [ -z "$1" ]
  then
    echo "No argument supplied. Please provide folder path to be cleaned up as a command line argument."
    exit 1	
elif [ -z "$2" ]
	then
    echo "Second argument not supplied. Please provide folder path to the trash folder as a command line argument."
    exit 1
fi

#change the working directory
cd $1

tmp_folder=$2
mkdir -p $tmp_folder

#clean up the results folder
for d in *; do
	mkdir $tmp_folder/$d
	mv $d/* $tmp_folder/$d/
	cp -R $tmp_folder/$d/outs $d
done
