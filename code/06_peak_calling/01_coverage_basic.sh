#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=100G
#SBATCH --job-name=01_coverage_basic
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL

## Explicitly pipe script output to a log
log_path=logs/01_coverage_basic.txt

{
set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module - Use this conda R version to keep Seurat-Object structure 
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_coverage_basic.R

echo "**** Job ends ****"
date

} > $log_path 2>&1


## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
