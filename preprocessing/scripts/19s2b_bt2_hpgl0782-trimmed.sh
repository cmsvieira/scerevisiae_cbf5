#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --workdir=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782
#SBATCH --partition=dpart
#SBATCH --qos=workstation
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=s2b_bt2_hpgl0782-trimmed
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/outputs/s2b_bt2_hpgl0782-trimmed.sbatchout

echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/scripts/19s2b_bt2_hpgl0782-trimmed.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782 || exit

## Converting the text sam to a compressed, sorted, indexed bamfile.
## Also printing alignment statistics to outputs/bowtie2_scerevisiae/hpgl0782.bam.stats
## This job depended on: 438907

if $(test ! -r outputs/bowtie2_scerevisiae/hpgl0782.sam); then
    echo "Could not find the samtools input file."
    exit 1
fi
samtools view -u -t /cbcbhomes/abelew/libraries/genome/scerevisiae.fasta \
  -S outputs/bowtie2_scerevisiae/hpgl0782.sam -o outputs/bowtie2_scerevisiae/hpgl0782.bam \
  2>outputs/bowtie2_scerevisiae/hpgl0782.bam.err 1>outputs/bowtie2_scerevisiae/hpgl0782.bam.out && \
  samtools sort -l 9 outputs/bowtie2_scerevisiae/hpgl0782.bam -o outputs/bowtie2_scerevisiae/hpgl0782-sorted.bam \
  2>outputs/bowtie2_scerevisiae/hpgl0782-sorted.err 1>outputs/bowtie2_scerevisiae/hpgl0782-sorted.out && \
  rm outputs/bowtie2_scerevisiae/hpgl0782.bam && \
  rm outputs/bowtie2_scerevisiae/hpgl0782.sam && \
  mv outputs/bowtie2_scerevisiae/hpgl0782-sorted.bam outputs/bowtie2_scerevisiae/hpgl0782.bam && \
  samtools index outputs/bowtie2_scerevisiae/hpgl0782.bam
## The following will fail if this is single-ended.
samtools view -b -f 2 -o outputs/bowtie2_scerevisiae/hpgl0782-paired.bam outputs/bowtie2_scerevisiae/hpgl0782.bam && \
  samtools index outputs/bowtie2_scerevisiae/hpgl0782-paired.bam
bamtools stats -in outputs/bowtie2_scerevisiae/hpgl0782.bam 2>outputs/bowtie2_scerevisiae/hpgl0782.bam.stats 1>&2 && \
  bamtools stats -in outputs/bowtie2_scerevisiae/hpgl0782-paired.bam 2>outputs/bowtie2_scerevisiae/hpgl0782-paired.stats 1>&2

## The following lines give status codes and some logging
echo $? > outputs/status/s2b_bt2_hpgl0782-trimmed.status
echo "###Finished ${SLURM_JOBID} 19s2b_bt2_hpgl0782-trimmed.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt

