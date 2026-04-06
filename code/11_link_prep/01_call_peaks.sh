#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=100G
#SBATCH --job-name=01_call_peaks
#SBATCH -c 1
#SBATCH -t 5-0:00:00
#SBATCH -o ../../processed-data/11_link_prep/logs/01_call_peaks.txt
#SBATCH -e ../../processed-data/11_link_prep/logs/01_call_peaks.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5

## List current modules for reproducibility
module list

Rscript 01_call_peaks.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/
