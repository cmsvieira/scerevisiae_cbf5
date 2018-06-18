#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --workdir=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782
#SBATCH --partition=dpart
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --time=18:00:00
#SBATCH --job-name=xzun_bt2_hpgl0782-trimmed
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/outputs/xzun_bt2_hpgl0782-trimmed.sbatchout

echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/scripts/16xzun_bt2_hpgl0782-trimmed.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782 || exit

## Compressing the sequences which failed to align against /cbcbhomes/abelew/libraries/genome/indexes/scerevisiae using options  --very-sensitive -L 14 

nice -n 20 xz -f -9e  /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/hpgl0782_forward-trimmed.fastq    2> /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/outputs/xzun_bt2_hpgl0782-trimmed.err \
   1> /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/outputs/xzun_bt2_hpgl0782-trimmed.out
## The following lines give status codes and some logging
echo $? > outputs/status/xzun_bt2_hpgl0782-trimmed.status
echo "###Finished ${SLURM_JOBID} 16xzun_bt2_hpgl0782-trimmed.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt

