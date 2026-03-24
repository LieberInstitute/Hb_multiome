#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=14_metacell_link_peaks
#SBATCH -c 1
#SBATCH -t 10-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-324%10

## Define loops and appropriately subset each variable for the array task ID
all_target_cell_type=(Astrocyte Endo Excit.Thal Inhib.Thal LHb.1 LHb.1.3 LHb.1.3.4 LHb.2.7 LHb.4 LHb.7 MHb.1 MHb.1.2 MHb.2 MHb.3 Microglia Oligo OPC Thal)
target_cell_type=${all_target_cell_type[$(( $SLURM_ARRAY_TASK_ID / 18 % 18 ))]}

all_other_cell_type=(Astrocyte Endo Excit.Thal Inhib.Thal LHb.1 LHb.1.3 LHb.1.3.4 LHb.2.7 LHb.4 LHb.7 MHb.1 MHb.1.2 MHb.2 MHb.3 Microglia Oligo OPC Thal)
other_cell_type=${all_other_cell_type[$(( $SLURM_ARRAY_TASK_ID / 1 % 18 ))]}

## Explicitly pipe script output to a log
log_path=../../processed-data/12_new_peaks/logs/14_metacell_link_peaks_${target_cell_type}_${other_cell_type}_${SLURM_ARRAY_TASK_ID}.txt

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

## Load the R module
module load conda_R/4.5

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 14_metacell_link_peaks.R --target_cell_type ${target_cell_type} --other_cell_type ${other_cell_type}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/

