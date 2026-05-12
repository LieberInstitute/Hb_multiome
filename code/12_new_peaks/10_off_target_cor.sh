#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=40G
#SBATCH --job-name=10_off_target_cor
#SBATCH -c 4
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/12_new_peaks/logs/10_off_target_cor.txt
#SBATCH -e ../../processed-data/12_new_peaks/logs/10_off_target_cor.txt

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

Rscript 10_off_target_cor.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/
