#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --workdir=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782
#SBATCH --partition=dpart
#SBATCH --qos=workstation
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --job-name=bcf_hpgl0782_scerevisiae
#SBATCH --mem=6G
#SBATCH --cpus-per-task=4
#SBATCH --output=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/outputs/bcf_hpgl0782_scerevisiae.sbatchout

echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/scripts/80bcf_hpgl0782_scerevisiae.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782 || exit

## Use samtools, bcftools, and vcfutils to get some idea about how many variant positions are in the data.
mkdir -p outputs/vcfutils_scerevisiae
echo "Started samtools sort at $(date)" >> outputs/vcfutils_scerevisiae.out
samtools sort -l 9 -@ 4 outputs/bowtie2_scerevisiae/hpgl0782.bam -o outputs/vcfutils_scerevisiae/hpgl0782.bam 2>outputs/samtools_sort.out 1>&2
if [ "$?" -ne "0" ]; then
    echo "samtools sort failed."
exit 1
fi

if [ ! -r "/cbcbhomes/abelew/libraries/genome/scerevisiae.fasta.fai" ]; then
    samtools faidx /cbcbhomes/abelew/libraries/genome/scerevisiae.fasta
fi
samtools mpileup -uvf /cbcbhomes/abelew/libraries/genome/scerevisiae.fasta 2>samtools_mpileup.err \
    outputs/vcfutils_scerevisiae/hpgl0782.bam |\
  bcftools call -c - 2>bcftools_call.err |\
  bcftools view -l 9 -o outputs/vcfutils_scerevisiae/hpgl0782.bcf -O b - \
    2>outputs/vcfutils_scerevisiae/hpgl0782_summary.err
if [ "$?" -ne "0" ]; then
    echo "mpileup/bcftools failed."
    exit 1
fi
bcftools index outputs/vcfutils_scerevisiae/hpgl0782.bcf 2>bcftools_index.err
echo "Successfully finished." >> outputs/vcfutils_scerevisiae.out

## The following lines give status codes and some logging
echo $? > outputs/status/bcf_hpgl0782_scerevisiae.status
echo "###Finished ${SLURM_JOBID} 80bcf_hpgl0782_scerevisiae.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt

