#!/bin/sh

##Submission script for the dada2 tutorial with run parameters that should send the job to c3
#PBS -m bea
#PBS -M btupper@bigelow.org
#PBS -N dada2-script
#PBS -q route

#PBS -l ncpus=24,mem=30gb
#PBS -l walltime=12:00:00

#PBS -e /home/btupper/edna/pbs_output
#PBS -o /home/btupper/edna/pbs_output

# Load modules
module use /mod/bigelow
module load dada2

path="/home/btupper/edna/edna-dada2"
version="v0.001"

cd ${path}

Rscript ${opts} ${path}/R/run_dada2_${version}.R ${path}/config/run_dada2_${version}.yml
