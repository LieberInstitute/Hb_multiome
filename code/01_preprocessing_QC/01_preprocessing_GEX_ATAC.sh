#!/bin/bash -l
#SBATCH --partition=katun	
#SBATCH --job-name=CR_QC_ini
#SBATCH --mem=50GB						                            # each job from the array will get its own private 30G to work with)
# SBATCH --mail-user=bioinformatic2019@gmail.com
# SBATCH --mail-type=end,fail
#SBATCH --array=1-3                                       #change the number 2 to the number of entries in array_targets.txt
id=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../array_target_names_2024a.txt)
#SBATCH --output=/logs/CR_QC_ini_%j.out	                        
#SBATCH --error=/logs/CR_QC_ini_%j.err	



echo "**** Job starts ****"
echo "     Build Seurat objects and calculate basic statistics for GEX and ATAC"
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

Rscript 01_preprocessing_GEX_ATAC.R $id

echo "== End of Job =="

## This script was made using JHPCE 3.0 SLURM Cluster
## Available from http://xxxx
## CSC. Jan 22nd, 2024 / 
