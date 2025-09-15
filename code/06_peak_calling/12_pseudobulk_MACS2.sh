#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=12_pseudobulk_MACS2
#SBATCH -c 4
#SBATCH -t 3-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL
#SBATCH --array=0-1%2   # 2 runs

set -eo pipefail

## Compute pseubobulk on multiome-seurat: regular and merged peaks
peak_ds=(regular merged)   # index 0..1
i=${SLURM_ARRAY_TASK_ID}
m=${#peak_ds[@]}

# guard
if (( i < 0 || i >= m )); then
  echo "Invalid task index: $i (must be 0..$((m-1)))"
  exit 1
fi

res="${peak_ds[$i]}"

# skip Broad
if [[ "$res" == "regular" ]]; then
  echo "[$(date)] Skipping peak_ds='regular' for array task $i"
  exit 0   # success so SLURM won’t retry
fi

#mkdir -p logs
log_path="logs/12_pseudobulk_MACS2_${res}_task${i}.log"

{

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module - compatible with seurat version + R_conda env
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 12_pseudobulk_MACS2.R --peak_ds "${res}"
ret=$?

echo "**** Job ends ****"
date
echo "Exit code: $ret"
exit $ret

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
