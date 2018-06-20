#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --workdir=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782
#SBATCH --partition=dpart
#SBATCH --qos=workstation
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=bt2_scerevisiae_hpgl0782-trimmed
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --output=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/outputs/bt2_scerevisiae_hpgl0782-trimmed.sbatchout

echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/scripts/15bt2_scerevisiae_hpgl0782-trimmed.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782 || exit

## This is a bowtie2 alignment of  /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/hpgl0782_forward-trimmed.fastq  against
## /cbcbhomes/abelew/libraries/genome/indexes/scerevisiae using arguments:  --very-sensitive -L 14 .
## This jobs depended on: 

if [[ -e "/scratch0" ]]; then
  scratchdir=$(mktemp -d "/scratch0/${USER}.XXXX")
  echo "Working in: ${scratchdir} on $(hostname)."
  cd "${scratchdir}" || exit
fi

mkdir -p outputs/bowtie2_scerevisiae && \
  sleep 3 && \
  bowtie2 -x /cbcbhomes/abelew/libraries/genome/indexes/scerevisiae  --very-sensitive -L 14  \
    -p 4 \
    -q   /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/hpgl0782_forward-trimmed.fastq  \
    --un outputs/bowtie2_scerevisiae/hpgl0782_unaligned_scerevisiae.fastq \
    --al outputs/bowtie2_scerevisiae/hpgl0782_aligned_scerevisiae.fastq \
    -S outputs/bowtie2_scerevisiae/hpgl0782.sam \
    2>outputs/bowtie2_scerevisiae/hpgl0782.err \
    1>outputs/bowtie2_scerevisiae/hpgl0782.out

if [[ -e "/scratch0" ]]; then
  cd /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782 && \
    rsync -a "${scratchdir}/outputs/bowtie2_scerevisiae/" outputs/bowtie2_scerevisiae && \
    rm -r "${scratchdir}"
fi

## The following lines give status codes and some logging
echo $? > outputs/status/bt2_scerevisiae_hpgl0782-trimmed.status
echo "###Finished ${SLURM_JOBID} 15bt2_scerevisiae_hpgl0782-trimmed.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt

