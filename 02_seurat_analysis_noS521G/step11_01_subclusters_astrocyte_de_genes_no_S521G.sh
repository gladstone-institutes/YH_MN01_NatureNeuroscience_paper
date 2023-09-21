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
mkdir -p /gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/results/data/11_subclusters_de_genes_no_S521G/astrocyte/

data_dir=/gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/
script_dir=/gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/scripts/GB-MN-1318/02_seurat_analysis_noS521G
container_dir=/gladstone/bioinformatics/projects/mn-1318-maxine-nelson-yadong-huang-snrnaseq-mm10-aug-2022/assets
export SINGULARITY_BINDPATH="$data_dir"

singularity exec $container_dir/gb_mn_1318_r_seurat4_harmony.sif Rscript $script_dir/src/11_01_subclusters_astrocyte_de_genes_no_S521G.R


## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
