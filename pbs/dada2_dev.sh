#!/bin/sh

#PBS -N dada2-dev
#PBS -q devel

#PBS -l ncpus=8,mem=8gb
#PBS -l walltime=8:00:00

#PBS -e /home/btupper/edna/pbs_output
#PBS -o /home/btupper/edna/pbs_output

# Load modules
module use /mod/bigelow
module load dada2

cd /home/btupper/edna/edna-dada2