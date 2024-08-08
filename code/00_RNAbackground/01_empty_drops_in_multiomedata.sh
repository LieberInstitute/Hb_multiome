#!/bin/bash
#SBATCH --partition=katun	
#SBATCH --job-name=empty_droplets
#SBATCH --mem=30GB						                                          # each job from the array will get its own private 30G to work with)
#SBATCH --array=1-5                                           	        #change the number 2 to the number of entries in array_targets.txt
id=$( sed -n ${SLURM_ARRAY_TASK_ID}p ../array_targets_names.txt)				                                      # partition or queue name
#SBATCH --output=/logs/ED_%j.out						                              
#SBATCH --error=/logs/ED_%j.err						                              

# NOTE: This script runs IN THE SAME DIRECTORY in which you ran sbatch
#       So include a cd command to ensure that you run it in the expected
#	directory (where your data files are located) OR use absolute paths
#	when specifying your input and output files

echo "**** Job starts ****"
echo "This script runs empty droplets in cellranger-ARC datasets"
echo "Processing cellranger-ARC experiments"
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
Rscript 01_empty_drops_in_multiomedata.R $id
echo "== End of Job =="

## Script for SLURM
## CSC. Jan 18th, 2024
## Last modif. Aug 08, 2024

## run command
## $ sbatch slurm_emptydrops_V2.sh

