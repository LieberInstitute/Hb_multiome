#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=02_GO_enrichment_cellType_universe
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-3%3

# Define an array mapping task ID to resolution level
# Task 1 -> broad
# Task 2 -> semi_broad
# Task 3 -> mid
RESOLUTIONS=( "broad" "semi_broad" "mid" )

# Select the resolution level for this specific task
# SLURM_ARRAY_TASK_ID is 1-based, so we use ${SLURM_ARRAY_TASK_ID}-1 for the 0-based array index
RES_LEVEL=${RESOLUTIONS[${SLURM_ARRAY_TASK_ID}-1]}

# Explicitly pipe script output to a log, using the resolution level in the name
log_path=logs/02_GO_enrichment_cellType_universe_${RES_LEVEL}.txt

{
set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "Resolution Level: $RES_LEVEL"

## Load the R module
module load conda_R/4.4.x

## List current modules for reproducibility
module list

Rscript 02_GO_enrichment_cellType_universe.R -r $RES_LEVEL

echo "**** Job ends ****"
date
} > $log_path 2>&1
