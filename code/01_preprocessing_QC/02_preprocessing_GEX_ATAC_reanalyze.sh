#!/bin/bash -l
#SBATCH --partition=katun	
#SBATCH --job-name=02_preprocessing_GEX_ATAC_reanalyze
#SBATCH --mem=50GB						                            # each job from the array will get its own private 30G to work with)
# SBATCH --mail-user=bioinformatic2019@gmail.com
# SBATCH --mail-type=end,fail
#SBATCH --array=1-10                                       #change the number 2 to the number of entries in array_targets.txt
#SBATCH --output=logs/02_preprocessing_GEX_ATAC_reanalyze_%A_%a.out
#SBATCH --error=logs/02_preprocessing_GEX_ATAC_reanalyze_%A_%a.err
id=$(sed -n ${SLURM_ARRAY_TASK_ID}p target_names_cellrangerARC_reanalyze.txt)


echo "**** Job starts ****"
echo "     Build Seurat objects and calculate basic statistics for GEX and ATAC from Cell Ranger Reanalyze"
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
module load conda_R/4.3.x
module list

echo "== This is the script =="

Rscript 02_preprocessing_GEX_ATAC_reanalyze.R $id

echo "== End of Job =="

## This script was made using JHPCE 3.0 SLURM Cluster
## Available from http://xxxx
## CSC. Jan 22nd, 2024 / 
