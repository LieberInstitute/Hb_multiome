#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=02_marker_compare
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/17_species_diverg/02_marker_compare/logs/02_marker_compare.txt
#SBATCH -e ../../processed-data/17_species_diverg/02_marker_compare/logs/02_marker_compare.txt

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

Rscript 02_marker_compare.R

echo "**** Job ends ****"
date
