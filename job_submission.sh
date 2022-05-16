#!/usr/bin/bash
#PBS -N run_dada2_pipeline
#PBS -q devel
#PBS -l ncpus=4,mem=5gb
#PBS -l walltime=2:00:00
#PBS -M user@bigelow.org
#PBS -m abe

module load dada2
data_path=/mnt/storage/data/edna/dada/projects/
Rscript $data_path"/asv_18S_preprocess.R" $data_path"/asv_18S_preprocess.yaml"

#Rscript $data_path"/asv_18S_postprocess.R" $data_path"/asv_18S_postprocess.yaml"


