#!/bin/bash
#$ -cwd
#$ -o /gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/tmp/
#$ -e /gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/tmp
#$ -pe smp 1
#$ -l mem_free=50G
#$ -l scratch=50G
#$ -l h_rt=06:00:00
#$ -j yes

#make the results output directory
mkdir -p /gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/results/data/06_subcluster_microglia_clusters_17_19/
mkdir -p /gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/results/plot/06_subcluster_microglia_clusters_17_19/

data_dir=/gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/
script_dir=/gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/scripts/GB-MN-1318/02_seurat_analysis_noS521G
container_dir=/gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/assets
export SINGULARITY_BINDPATH="$data_dir"

singularity exec $container_dir/gb_mn_1318_r_seurat4_harmony.sif Rscript $script_dir/src/06_subcluster_microglia_clusters_17_19.R


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"

