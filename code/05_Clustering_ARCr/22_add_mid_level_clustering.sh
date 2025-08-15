#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=22_add_mid_level_clustering
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/off/null
#SBATCH -e /dev/off/null
# SBATCH --mail-type=ALL

set -eo pipefail

{

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.3.x # kep chromatin obj compatible

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 22_add_mid_level_clustering.R
# Capture return code and exit safely
ret=$?

echo "**** Job ends ****"
date
echo "Exit code: $ret"
exit $ret

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
