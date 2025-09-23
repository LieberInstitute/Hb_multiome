#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=30G
#SBATCH --job-name=18_pseudobulk_DARs_Volcano_Violin
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o logs/18_pseudobulk_DARs_Volcano_Violin.txt
#SBATCH -e logs/18_pseudobulk_DARs_Volcano_Violin.txt
#SBATCH --mail-type=ALL

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
module load conda_R/4.4

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 18_pseudobulk_DARs_Volcano_Violin.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
