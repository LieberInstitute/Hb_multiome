#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=01_hashikawa_markers
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/17_species_diverg/01_hashikawa_markers/logs/01_hashikawa_markers.txt
#SBATCH -e ../../processed-data/17_species_diverg/01_hashikawa_markers/logs/01_hashikawa_markers.txt

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

Rscript 01_hashikawa_markers.R

echo "**** Job ends ****"
date