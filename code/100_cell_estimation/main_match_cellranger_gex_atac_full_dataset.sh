#!/bin/bash
#SBATCH --partition=katun	
#SBATCH --job-name=main_match_cellranger_gex_atac_full_dataset
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=15GB						                                    # each job from the array will get its own private 30G to work with)
# SBATCH --array=1-12                                           	  # change the array range acordingly with the rows listed in the array_targets.txt
# id=$( sed -n ${SLURM_ARRAY_TASK_ID}p ../array_target_names_full_dataset.txt)				                                  
#SBATCH -o logs/main_match_cellranger_gex_atac_full_dataset.txt
#SBATCH -e logs/main_match_cellranger_gex_atac_full_dataset.txt
# SBATCH --mail-type=ALL

set -e

echo "**** Job starts ****"
echo "General script to run CellRanger-ARC cell matching"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job name (script): ${SLURM_JOB_NAME}"
echo "Hostname (computer node): ${HOSTNAME}"
echo ""
echo "Job id: ${SLURM_JOBID}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "Array/sample: $id"

## load modules
module conda_R/4.3.x
## List current modules for reproducibility
module list

MAINDIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome"
CODEDIR="${MAINDIR}/code"
PROCESSEDIR="${MAINDIR}/processed-data"
# PLOTDIR="${MAINDIR}/plots"


## Delete the logs/old-results, and re-submit seurat builder
cd ${CODEDIR}/100_cell_estimation
# mkdir -p logs ## Create the logs directory if it doesn't exist
# rm logs/*.txt
# mv ${PROCESSEDIR}/100_cell_match_cellranger_gex_atac/*.tsv tmp/
# mv ${PROCESSEDIR}/100_cell_match_cellranger_gex_atac/*.csv tmp/

## Main job-array
id1=$(sbatch --parsable 01_cell_match_cellranger_gex_atac.sh)

## Dependency job $id1
sbatch --dependency=afterok:$id1 02_cell_match_stats_integrated.sh

echo "Parsable-job ID dependency: ${id1}"
echo "Completed!"

echo "**** Job ends ****"

date




