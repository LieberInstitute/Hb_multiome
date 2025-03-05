#!/bin/bash
#SBATCH --partition=katun
#SBATCH --mem=30G
#SBATCH --job-name=02_integrate_seurats_CRarc
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=1-2%20

## Define loops and appropriately subset each variable for the array task ID
all_type_mtx=(data_counts norm_counts)
type_mtx=${all_type_mtx[$(( $SLURM_ARRAY_TASK_ID / 1 % 2 ))]}

## Explicitly pipe script output to a log
log_path=logs/02_integrate_seurats_CRarc_${type_mtx}_${SLURM_ARRAY_TASK_ID}.txt

{
set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "Assay type: ${type_mtx}"

## Load the R module
module load conda_R/4.3.x
module list

## Edit with your job command
Rscript 02_integrate_seurats_CRarc.R --type_mtx ${type_mtx}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/

