#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=30G
#SBATCH --job-name=02_Hb_celltypes_from_seurat_reanalyze
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o logs/02_Hb_celltypes_from_seurat_reanalyze.%a.txt
#SBATCH -e logs/02_Hb_celltypes_from_seurat_reanalyze.%a.txt
# SBATCH --mail-type=ALL
#SBATCH --array=1-8%20
id=$(sed -n ${SLURM_ARRAY_TASK_ID}p input_wnn_rds_names.txt)

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "WNN file name: $id"

## Load the R module
module load conda_R/4.4.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02_Hb_celltypes_from_seurat_reanalyze.R $id

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/
