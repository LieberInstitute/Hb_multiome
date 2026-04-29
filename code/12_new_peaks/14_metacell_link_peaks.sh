#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=14_metacell_link_peaks
#SBATCH -c 1
#SBATCH -t 10-0:00:00
#SBATCH -o ../../processed-data/12_new_peaks/logs/14_metacell_link_peaks_%a.txt
#SBATCH -e ../../processed-data/12_new_peaks/logs/14_metacell_link_peaks_%a.txt
#SBATCH --array=1-17%17

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

Rscript 14_metacell_link_peaks.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/
