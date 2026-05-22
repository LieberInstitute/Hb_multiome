#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=150G
#SBATCH --job-name=10_celltype_func_annot
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/05_03_annotation_adjustments/10_celltype_func_annot/logs/10_celltype_func_annot.txt
#SBATCH -e ../../processed-data/05_03_annotation_adjustments/10_celltype_func_annot/logs/10_celltype_func_annot.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"

## Load the R module
module load conda_R/4.5

## List current modules for reproducibility
module list

Rscript 10_celltype_func_annot.R

echo "**** Job ends ****"
date
