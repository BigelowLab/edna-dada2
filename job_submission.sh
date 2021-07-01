#!/usr/bin/bash
#PBS -N run_dada2_pipeline
#PBS -q normal
#PBS -l ncpus=16,mem=50gb
#PBS -l walltime=72:00:00
#PBS -M user@bigelow.org
#PBS -m abe

module load dada2
data_path=/mnt/storage/data/edna/dada/projects/
Rscript $data_path"/dada2_18S.R" $data_path"/dada2_18S.yaml"

