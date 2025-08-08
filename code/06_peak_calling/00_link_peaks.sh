#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=00_link_peaks
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL
#SBATCH --array=0-5%6   # 2 methods * 3 windows = 6 tasks

# 0: pearson, 1e5
# 1: pearson, 5e4
# 2: spearman, 1e5
# 3: spearman, 5e4

set -eo pipefail

peak_methods=(pearson spearman)   # m=2
window_size=(1e5 5e4 2.5e4)       # n=3 (added 2.5e4)

i=${SLURM_ARRAY_TASK_ID}
m=${#peak_methods[@]}
n=${#window_size[@]}
total=$(( m * n ))

if (( i < 0 || i >= total )); then
  echo "Invalid task index: $i (total=$total)"; exit 1
fi

method_idx=$(( i / n ))
win_idx=$(( i % n ))

p_met="${peak_methods[$method_idx]}"
w_size="${window_size[$win_idx]}"

#mkdir -p logs
log_path="logs/00_link_peaks_method_${p_met}_window_${w_size}_task_${i}.log"

{

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 00_link_peaks.R --p_met "${p_met}" --w_size "${w_size}"

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
