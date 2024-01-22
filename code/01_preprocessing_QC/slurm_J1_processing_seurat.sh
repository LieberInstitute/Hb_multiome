#!/bin/bash -l
#SBATCH --mail-user=bioinformatic2019@gmail.com
#SBATCH --mail-type=end,fail
#SBATCH --output=logs/cellranger01_J1_%A_%a.out                             # Standard Output File. %A is replaced by SLURM_ARRAY_JOB_ID and %a by SLURM_ARRAY_TASK_ID 
#SBATCH --error=logs/cellranger01_J1_%A_%a.err   
#SBATCH --array=1-2                                         	    #change the number 2 to the number of entries in array_targets.txt
id=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../array_targets_names.txt)
#SBATCH --partition=shared					                        # partition or queue name
#SBATCH --job-name=cellranger0   			                        # name job for easier spotting, controlling
#SBATCH --mem=50GB						                            # each job from the array will get its own private 30G to work with)


echo "**** Job starts ****"
echo "     Build Seurat objects and calculate statistics for GEX and ATAC"
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
module load R
# list modules
module list

echo "== This is the script =="

Rscript 01_preprocessing_GEX_ATAC.R $id

echo "== End of Job =="

## This script was made using JHPCE 3.0 SLURM Cluster
## Available from http://xxxx
## CSC. Jan 22nd, 2024 / 
