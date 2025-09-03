#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=100G
#SBATCH --job-name=02_link_peaks_MACS2
#SBATCH -c 2
#SBATCH -t 2-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL
#SBATCH --array=0-1%2   # 2 methods

set -eo pipefail

clust_res=(Broad Mid)   # index 0..1

i=${SLURM_ARRAY_TASK_ID}
m=${#clust_res[@]}

# guard
if (( i < 0 || i >= m )); then
  echo "Invalid task index: $i (must be 0..$((m-1)))"
  exit 1
fi

res="${clust_res[$i]}"

# # skip Broad
# if [[ "$res" == "Broad" ]]; then
#   echo "[$(date)] Skipping clust_res='Broad' for array task $i"
#   exit 0   # success so SLURM won’t retry
# fi

#mkdir -p logs
log_path="logs/02_link_peaks_MACS2_${res}_task_${i}.log"

{

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "Selected clust_res: ${res}"

## Load the R module
module load conda_R/4.3.x

## List current modules for reproducibility
module list
which Rscript || true

## Edit with your job command
Rscript 02_link_peaks_MACS2.R --clust_res "${res}"
ret=$?

echo "**** Job ends ****"
date
echo "Exit code: $ret"
exit $ret

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/