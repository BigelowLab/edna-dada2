##Submission script for the dada2 tutorial with run parameters that should send the job to c3
#PBS -m bea
#PBS -M btupper@bigelow.org
#PBS -N dada2-script
#PBS -q route

#PBS -l ncpus=24,mem=30gb
#PBS -l walltime=12:00:00

#PBS -e /home/pcountway/scripts/out/SK18S_batch
#PBS -o /home/pcountway/scripts/out/SK18S_batch

# Load modules
module use /mod/bigelow
module load dada2
path=/home/btupper/edna/edna-dada2
opts="--no-save --no-restore --no-site-file --no-environ"
Rscript ${opts} ${path}/Rscript/run_dada2_v0.000.R ${path}/config/run_dada2_v0.000.yml
