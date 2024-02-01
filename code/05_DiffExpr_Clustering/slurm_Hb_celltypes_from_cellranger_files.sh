#!/bin/bash -l
#SBATCH --mail-user=bioinformatic2019@gmail.com
#SBATCH --mail-type=end,fail
#SBATCH --output=logs/S05_%A_%a.out                             # Standard Output File. %A is replaced by SLURM_ARRAY_JOB_ID and %a by SLURM_ARRAY_TASK_ID 
#SBATCH --error=logs/S05_%A_%a.err   
#SBATCH --job-name=S05_celltype_arc		
#SBATCH --partition=shared					                        # partition or queue name
#SBATCH --mem=20G                                                  # ram per node	
#SBATCH --cpus-per-task=2                                           # number of CPUs per process					                           
#SBATCH --array=1-2                                         	    #change the number 2 to the number of entries in array_targets.txt
id=$(sed -n ${SLURM_ARRAY_TASK_ID}p ../array_targets_names.txt)

echo "**** Job starts ****"
echo "     Run Gene Markers Identification for cellranger-ARC datasets"
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

Rscript 05_Hb_celltypes_from_cellranger_files.R $id

echo "== End of Job =="

## This script was made using JHPCE 3.0 SLURM Cluster
## Available from http://xxxx
## CSC. Feb 1st, 2024
