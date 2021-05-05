#!/usr/bin/bash
#PBS -N run_dada2_pipeline
#PBS -q normal
#PBS -l ncpus=16,mem=50gb
#PBS -l walltime=72:00:00
#PBS -M rsleith@bigelow.org
#PBS -m abe

module load dada2
data_path=/mnt/storage/data/edna/dada/projects/robin_foo/may
#Rscript /mnt/storage/data/edna/dada/projects/robin_foo/highland/dada.R /mnt/storage/data/edna/dada/projects/robin_foo/highland/dada18S.yaml
Rscript $data_path"/dada18S.R" $data_path"/dada18S.yaml"

