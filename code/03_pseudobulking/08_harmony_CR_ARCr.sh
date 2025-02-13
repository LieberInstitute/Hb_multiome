#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=60G
#SBATCH --job-name=08_harmony_CR_ARCr
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL

## Explicitly pipe script output to a log
log_path=logs/08_harmony_CR_ARCr_${SLURM_JOB_ID}.txt

{
set -e

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
Rscript 08_harmony_CR_ARCr.R

echo "**** Job ends ****"
date

} > $log_path 2>&1


## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
