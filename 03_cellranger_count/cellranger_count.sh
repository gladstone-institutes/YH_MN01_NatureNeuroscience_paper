#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/tmp/
#$ -e /gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/tmp/
#$ -pe smp 4
#$ -l mem_free=40G
#$ -l scratch=50G
#$ -l h_rt=08:00:00
#$ -j yes
#$ -N cellranger
#$ -V

#set the paths
base_dir=/gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022
transcriptome_dir=/gladstone/bioinformatics/adnetworksppg/Project_1_Huang/NK01/assets/reference_genomes/analysis_hapoe_chr

#change the working directory
cd $base_dir/results/cellranger_count/

#process command line arguments to the script
sample_id=$1
sample_name=$2

#load local modules
module load CBI
module load cellranger/7.0.0

#run cell ranger count
cellranger count \
--id=$sample_id \
--transcriptome=$transcriptome_dir/adppg-mm10-apoe-chr-mapt-chr \
--fastqs=$base_dir/data/fastq/YH_MN02 \
--sample=$sample_name \
--include-introns=true \
--localcores=32


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
