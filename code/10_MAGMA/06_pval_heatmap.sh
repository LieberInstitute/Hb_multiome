#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=06_pval_heatmap
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-33%10

## Define loops and appropriately subset each variable for the array task ID
all_gwas=(MDD MDD2019 panic SCZ SUD2020 AUD CUD ext_cannabis lifetime_cannabis OUD SUD2)
gwas=${all_gwas[$(( $SLURM_ARRAY_TASK_ID / 3 % 11 ))]}

all_cell_type_group=(broad semi_broad mid)
cell_type_group=${all_cell_type_group[$(( $SLURM_ARRAY_TASK_ID / 1 % 3 ))]}

## Explicitly pipe script output to a log
log_path=../../processed-data/10_MAGMA/logs/06_pval_heatmap_${gwas}_${cell_type_group}_${SLURM_ARRAY_TASK_ID}.txt

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
Rscript 06_pval_heatmap.R --gwas ${gwas} --cell_type_group ${cell_type_group}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/
