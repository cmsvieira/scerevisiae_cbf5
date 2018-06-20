#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --workdir=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782
#SBATCH --partition=dpart
#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=htall_bt2_hpgl0782-trimmed
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/outputs/htall_bt2_hpgl0782-trimmed.sbatchout

echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/scripts/21htall_bt2_hpgl0782-trimmed.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782 || exit

## Counting the number of hits in outputs/bowtie2_scerevisiae/hpgl0782.bam for each feature found in /cbcbhomes/abelew/libraries/genome/scerevisiae.gff
## Is this stranded? no.  The defaults of htseq are:
##  --order=name --idattr=gene_id --minaqual=10 --type=exon --stranded=yes --mode=union 


htseq-count  --help 2>&1 | tail -n 3
htseq-count -q -f bam -s no  --idattr exon_number  --type exon \
  outputs/bowtie2_scerevisiae/hpgl0782.bam \
  /cbcbhomes/abelew/libraries/genome/scerevisiae.gff \
  2>outputs/bowtie2_scerevisiae/hpgl0782_htseq.err \
  1>outputs/bowtie2_scerevisiae/hpgl0782.count && \
    xz -f -9e outputs/bowtie2_scerevisiae/hpgl0782.count 2>outputs/bowtie2_scerevisiae/hpgl0782_htseq.err.xz 1>outputs/bowtie2_scerevisiae/hpgl0782.count.xz

## The following lines give status codes and some logging
echo $? > outputs/status/htall_bt2_hpgl0782-trimmed.status
echo "###Finished ${SLURM_JOBID} 21htall_bt2_hpgl0782-trimmed.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt

