#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --mail-type=NONE
#SBATCH --workdir=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782
#SBATCH --partition=dpart
#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --job-name=bt2st_bt2_hpgl0782-trimmed
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/outputs/bt2st_bt2_hpgl0782-trimmed.sbatchout

echo "####Started /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/scripts/18bt2st_bt2_hpgl0782-trimmed.sh at $(date) on $(hostname)" >> outputs/log.txt
cd /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782 || exit

## This is a stupidly simple job to collect alignment statistics.

if [ ! -e "outputs/bowtie2_stats.csv" ]; then
    echo "original reads, single hits, failed reads, multi-hits, rpm" > outputs/bowtie2_stats.csv
fi
original_reads_tmp=$(grep " reads; of these" "outputs/bowtie2_scerevisiae/hpgl0782.err" 2>/dev/null | awk '{print $1}' | sed 's/ //g')
original_reads=${original_reads_tmp:-0}
one_align_tmp=$(grep " aligned exactly 1 time" "outputs/bowtie2_scerevisiae/hpgl0782.err" | awk '{print $1}' | sed 's/ .*//g')
one_align=${one_align_tmp:-0}
failed_tmp=$(grep " aligned 0 times" "outputs/bowtie2_scerevisiae/hpgl0782.err" | awk '{print $1}' | sed 's/ .*//g')
failed=${failed_tmp:-0}
sampled_tmp=$(grep " aligned >1 times" "outputs/bowtie2_scerevisiae/hpgl0782.err" | awk '{print $1}' | sed 's/ .*//g')
sampled=${sampled_tmp:-0}
rpm_tmp=$(perl -e "printf(1000000 / $(( ${one_align} + ${sampled} )) ) " 2>/dev/null)
rpm=${rpm_tmp:-0}
stat_string=$(printf "hpgl0782,v0M1,%s,%s,%s,%s,%s" "${original_reads}" "${one_align}" "${failed}" "${sampled}" "${rpm}")
echo "$stat_string" >> outputs/bowtie2_stats.csv
## The following lines give status codes and some logging
echo $? > outputs/status/bt2st_bt2_hpgl0782-trimmed.status
echo "###Finished ${SLURM_JOBID} 18bt2st_bt2_hpgl0782-trimmed.sh at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt

walltime=$(scontrol show job "${SLURM_JOBID}" | grep RunTime | perl -F'/\s+|=/' -lane '{print $F[2]}')
echo "#### walltime used by ${SLURM_JOBID} was: ${walltime:-null}" >> outputs/log.txt
maxmem=$(sstat --format=MaxVMSize -n "${SLURM_JOBID}.batch")
echo "#### maximum memory used by ${SLURM_JOBID} was: ${maxmem:-null}" >> outputs/log.txt
avecpu=$(sstat --format=AveCPU -n "${SLURM_JOBID}.batch")
echo "#### average cpu used by ${SLURM_JOBID} was: ${avecpu:-null}" >> outputs/log.txt

