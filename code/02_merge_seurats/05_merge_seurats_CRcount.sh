#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=100G
#SBATCH --job-name=05_merge_seurats_CRcount
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL
#SBATCH --array=1-2%20

## Define loops and appropriately subset each variable for the array task ID
all_type_mtx=(data_counts norm_counts)
type_mtx=${all_type_mtx[$(( $SLURM_ARRAY_TASK_ID / 1 % 2 ))]}

## Explicitly pipe script output to a log
log_path=logs/05_merge_seurats_CRcount_${type_mtx}_${SLURM_ARRAY_TASK_ID}.txt

{
set -e

echo "**** Job starts ****"
date

# echo "Removing previous Dir/Subdir/files outputs if exists"
# echo " "
# SUBDIR="02_merge_seurats/cellranger_count"
# MAINDIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome"
# PROCESSEDIR="${MAINDIR}/processed-data/${SUBDIR}"
# PROCESSEDIR_CSV="${MAINDIR}/processed-data/${SUBDIR}/csv_files"
# PLOTDIR="${MAINDIR}/plots/${SUBDIR}"

## Delete the logs/old-results, and re-submit seurat builder
# mkdir -p logs ## Create the logs directory if it doesn't exist
# rm logs/05_merge_seurats_CRcount_*.txt
# rm ${PROCESSEDIR}/seurat.*.rds
# rm ${PROCESSEDIR_CSV}/*.csv
# rm ${PLOTDIR}/*.png
echo " "

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.3.x
module list

## Edit with your job command
Rscript 05_merge_seurats_CRcount.R --type_mtx ${type_mtx}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/

