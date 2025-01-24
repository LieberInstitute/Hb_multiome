#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=02_v2_Hb_celltypes_ARCr
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL
#SBATCH --array=1-8%20
id=$(sed -n ${SLURM_ARRAY_TASK_ID}p input_wnn_rds_names.txt)

## Explicitly pipe script output to a log
log_path=logs/02_v2_Hb_celltypes_ARCr_${id}.txt

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
echo "WNN file name: $id"

## Load the R module
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02_Hb_celltypes_from_seurat_reanalyze_v2.R $id

echo "**** Job ends ****"
date

} > $log_path 2>&1


