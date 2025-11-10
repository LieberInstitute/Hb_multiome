#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=40G
#SBATCH --job-name=09_GO_Enrich
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-3%3  # Run 3 tasks (1=broad, 2=semi_broad, 3=mid)

# Define an array mapping task ID to resolution level
# Task 1 -> broad
# Task 2 -> semi_broad
# Task 3 -> mid
RESOLUTIONS=( "broad" "semi_broad" "mid" )

# Select the resolution level for this specific task
# SLURM_ARRAY_TASK_ID is 1-based, so we use ${SLURM_ARRAY_TASK_ID}-1 for the 0-based array index
RES_LEVEL=${RESOLUTIONS[${SLURM_ARRAY_TASK_ID}-1]}

# Explicitly pipe script output to a log, using the resolution level in the name
log_path=logs/09_GO_enrich_${RES_LEVEL}.txt

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

# Execute the R script, passing the resolution level using the -r argument
# NOTE: Replace 'your_go_enrichment_script.R' with your actual R file name
Rscript your_go_enrichment_script.R -r $RES_LEVEL

echo "**** Job ends ****"
date
} > $log_path 2>&1
