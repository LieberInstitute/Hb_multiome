#!/bin/bash
#SBATCH --partition=katun	
#SBATCH --job-name=match_cells
#SBATCH --mem=10GB						                                          # each job from the array will get its own private 30G to work with)
#SBATCH --array=1-10                                           	        #change the number 2 to the number of entries in array_targets.txt
id=$( sed -n ${SLURM_ARRAY_TASK_ID}p array_target_names.txt)   
#SBATCH --output=/logs/match_%j.out						                              
#SBATCH --error=/logs/match_%j.err						                              

# NOTE: This script runs IN THE SAME DIRECTORY in which you ran sbatch
#       So include a cd command to ensure that you run it in the expected
#	directory (where your data files are located) OR use absolute paths
#	when specifying your input and output files

echo "**** Job starts ****"
echo "This script gets the matching barcodes of multi-ome cell estimation between cellranger-count and cellranger-atac"
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
module list

echo "== This is the script =="
Rscript 01_cell_match_cellranger_gex_atac.R $id
echo "== End of Job =="

## Script for SLURM
## CSC. Sep, 2024
## run command
## $ sbatch file_name

