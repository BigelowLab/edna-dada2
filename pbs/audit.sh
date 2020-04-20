#PBS -m bea
#PBS -M btupper@bigelow.org
#PBS -N audit-script
#PBS -q route

#PBS -l select=1:ncpus=1
#PBS -l walltime=1:00:00

#PBS -e /home/btupper/pbs_out
#PBS -o /home/btupper/pbs_out

module use /mod/bigelow
module load dada2
path=/home/btupper/edna/edna-dada2
Rscript ${path}/Rscript/audit.R /home/btupper/edna/data/audits