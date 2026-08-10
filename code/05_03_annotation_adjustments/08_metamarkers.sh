#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=08_metamarkers
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o ../../processed-data/05_03_annotation_adjustments/logs/08_metamarkers.txt
#SBATCH -e ../../processed-data/05_03_annotation_adjustments/logs/08_metamarkers.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"

module load conda_R/4.5

## List current modules for reproducibility
module list

Rscript 08_metamarkers.R

echo "**** Job ends ****"
date
