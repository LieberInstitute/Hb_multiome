#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=02_cell_types_percentages_v3
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL
#SBATCH --array=1-8%20

## Define loops and appropriately subset each variable for the array task ID
all_cellranger_pipe=(CR_crossBarcodes CR_complementBarcodes CR_arc_reanalyze CR_arc_reanalyze_outliers)
cellranger_pipe=${all_cellranger_pipe[$(( $SLURM_ARRAY_TASK_ID / 2 % 4 ))]}

all_marker_lst=(literature_base data_driven)
marker_lst=${all_marker_lst[$(( $SLURM_ARRAY_TASK_ID / 1 % 2 ))]}

## Explicitly pipe script output to a log
log_path=logs/02_cell_types_percentages_v3_${cellranger_pipe}_${marker_lst}_${SLURM_ARRAY_TASK_ID}.txt

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

## Load the R module
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02_cell_types_percentages.R --cellranger_pipe ${cellranger_pipe} --marker_lst ${marker_lst}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/

