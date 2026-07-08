#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=25G
#SBATCH --job-name=01_halo_data
#SBATCH -c 6
#SBATCH -t 3-00:00:00
#SBATCH -o ../../processed-data/14_vgat_rnascope/01_halo_data/logs/01_halo_data.txt
#SBATCH -e ../../processed-data/14_vgat_rnascope/01_halo_data/logs/01_halo_data.txt

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

Rscript 01_halo_data.R

echo "**** Job ends ****"
date