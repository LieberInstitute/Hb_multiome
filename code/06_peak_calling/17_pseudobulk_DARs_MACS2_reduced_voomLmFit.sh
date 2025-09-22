#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=30G
#SBATCH --job-name=17_pseudobulk_DARs_MACS2_reduced_voomLmFit_ct
#SBATCH -c 4
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=0-17%20   # 18 clusters, max 4 running concurrently
# SBATCH --mail-type=ALL

set -eo pipefail

# Bash array - WNN at Mid level
clust=(
  "Astrocyte" "Endo" "Excit.Thal" "Inhib.Thal" "LHb.1" "LHb.1.3" "LHb.1.3.4"
  "LHb.2.7" "LHb.4" "LHb.7" "MHb.1" "MHb.1.2" "MHb.2" "MHb.3"
  "Microglia" "Oligo" "OPC" "Thal"
)

# Allow local testing; SLURM sets this in the array
i=${SLURM_ARRAY_TASK_ID:-0}
m=${#clust[@]}

# guard
if (( i < 0 || i >= m )); then
  echo "Invalid task index: $i (must be 0..$((m-1)))" >&2
  exit 1
fi

res="${clust[$i]}"

mkdir -p logs
log_path="logs/17_pseudobulk_DARs_MACS2_reduced_voomLmFit_ct${res}.log"

{

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 17_pseudobulk_DARs_MACS2_reduced_voomLmFit.R "${res}"
ret=$?

echo "**** Job ends ****"
date
echo "Exit code: $ret"
exit $ret

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
