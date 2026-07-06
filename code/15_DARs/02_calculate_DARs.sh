#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=8G
#SBATCH --job-name=02_calculate_DARs
#SBATCH -c 1
#SBATCH -t 10-00:00:00
#SBATCH -o ../../processed-data/15_DARs/logs/02_calculate_DARs_%a.txt
#SBATCH -e ../../processed-data/15_DARs/logs/02_calculate_DARs_%a.txt
#SBATCH --array=8-10%3
#SBATCH --exclude=compute-094

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
module load conda_R/4.5

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02_calculate_DARs.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.4.0
## available from http://research.libd.org/slurmjobs/
