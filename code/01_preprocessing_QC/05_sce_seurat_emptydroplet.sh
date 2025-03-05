#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=50G
#SBATCH --job-name=05_sce_seurat_emptydroplet
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o logs/05_sce_seurat_emptydroplet.%a.txt
#SBATCH -e logs/05_sce_seurat_emptydroplet.%a.txt
# SBATCH --mail-type=ALL
#SBATCH --array=1-10%20
id=$(sed -n ${SLURM_ARRAY_TASK_ID}p target_names_cellrangerARC_reanalyze.txt)

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo " "
echo "Sample: $id"

## Load the R module
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 05_sce_seurat_emptydroplet.R $id

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
