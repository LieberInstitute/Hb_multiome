#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=05_rename_idents
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --mail-type=ALL
#SBATCH --array=1-4%20

## Define loops and appropriately subset each variable for the array task ID

# Old rds objects with only rna harmonized and lsi
# all_clustering_name=("seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1" 
# "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvainM_lsi_r1" "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.leiden_lsi_r1" "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.SLM_lsi_r1")

# Old rds objects with only rna and lsi harmonized
all_clustering_name=("seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r1" 
"seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2" 
"seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r1" 
"seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r2")

clustering_name=${all_clustering_name[$(( $SLURM_ARRAY_TASK_ID / 1 % 4 ))]}

## Explicitly pipe script output to a log
log_path=logs/05_rename_idents_${clustering_name}_${SLURM_ARRAY_TASK_ID}.txt

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
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 05_rename_idents.R --clustering_name ${clustering_name}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/

