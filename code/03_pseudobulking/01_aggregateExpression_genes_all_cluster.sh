#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=01_aggregateExpression_genes_all_cluster
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL
#SBATCH --array=1-4%20

## Define loops and appropriately subset each variable for the array task ID
all_type_mtx=(data_counts norm_counts)
type_mtx=${all_type_mtx[$(( $SLURM_ARRAY_TASK_ID / 2 % 2 ))]}

all_integration_model=(CCA Harmony)
integration_model=${all_integration_model[$(( $SLURM_ARRAY_TASK_ID / 1 % 2 ))]}

## Explicitly pipe script output to a log
log_path=logs/01_aggregateExpression_genes_all_cluster_${type_mtx}_${integration_model}_${SLURM_ARRAY_TASK_ID}.txt

{
set -e

echo "**** Job starts ****"
date

echo "Removing previous Dir/Subdir/files outputs if exists"
echo " "
## remove previous outputs if exists
MAINDIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome"
PROCESSEDIR="${MAINDIR}/processed-data/03_pseudobulking/cellrangerARC_reanalyze"
PROCESSEDIR_CSV="${MAINDIR}/processed-data/03_pseudobulking/cellrangerARC_reanalyze/cvs_files_markers"
PLOTDIR="${MAINDIR}/plots/03_pseudobulking/cellrangerARC_reanalyze"

## Delete the logs/old-results, and re-submit seurat builder
# mkdir -p logs ## Create the logs directory if it doesn't exist
rm logs/01_aggregateExpression_genes_all_cluster*.txt
rm ${PROCESSEDIR}/seurat.*Harmony*.rds
rm ${PROCESSEDIR}/seurat.*CCA*.rds
rm ${PROCESSEDIR_CSV}/*.csv
rm ${PLOTDIR}/*.png
echo " "

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "Processing: " ${type_mtx}_${integration_model}

## Load the R module
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_aggregateExpression_genes_all_cluster.R --type_mtx ${type_mtx} --integration_model ${integration_model}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/

