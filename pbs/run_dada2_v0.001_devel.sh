#!/bin/sh

#PBS -N dada2-script
#PBS -q devel

#PBS -l ncpus=8,mem=8gb
#PBS -l walltime=8:00:00

#PBS -e /home/btupper/edna/pbs_output
#PBS -o /home/btupper/edna/pbs_output

# Load modules
module use /mod/bigelow
module load dada2

#path="/home/btupper/edna/edna-dada2"
#version="v0.001"

#cd ${path}

#Rscript ${opts} ${path}/R/run_dada2_${version}.R ${path}/config/run_dada2_${version}.yml
