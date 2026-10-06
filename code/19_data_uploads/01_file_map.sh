#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=01_file_map
#SBATCH -c 1
#SBATCH -t 20-00:00:00
#SBATCH -o ../../processed-data/19_data_uploads/logs/01_file_map.txt
#SBATCH -e ../../processed-data/19_data_uploads/logs/01_file_map.txt

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
module load conda_R/4.6

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_file_map.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.4.0
## available from http://research.libd.org/slurmjobs/
