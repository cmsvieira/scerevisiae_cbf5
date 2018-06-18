#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --workdir=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782
#SBATCH --partition=dpart
#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=parsesnp_scerevisiae
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/outputs/parsesnp_scerevisiae.sbatchout

echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/scripts/81parsesnp_scerevisiae.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782 || exit


## Parse the SNP data and generate a modified scerevisiae genome.
##  This should read the file:
## outputs/vcfutils_scerevisiae/hpgl0782.bcf
##  and provide 4 new files:
## outputs/vcfutils_scerevisiae/hpgl0782_scerevisiae_count.txt
## outputs/vcfutils_scerevisiae/hpgl0782_scerevisiae_all.txt
## outputs/vcfutils_scerevisiae/hpgl0782_scerevisiae_pct.txt
## and a modified genome: outputs/vcfutils_scerevisiae/hpgl0782_scerevisiae_modified.fasta

/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/scripts/81parsesnp_scerevisiae.pl

## The following lines give status codes and some logging
echo $? > outputs/status/parsesnp_scerevisiae.status
echo "###Finished ${SLURM_JOBID} 81parsesnp_scerevisiae.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt

