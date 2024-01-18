#!/bin/bash -l
#SBATCH --array=1-5                                           	        #change the number 2 to the number of entries in array_targets.txt
# sample_name=$( sed -n ${SLURM_ARRAY_TASK_ID}p array_targets.txt)      	# for each array job, I pull out the value from array_targets corresponding to that job
id=$( sed -n ${SLURM_ARRAY_TASK_ID}p array_targets_names.txt)
#SBATCH --partition=shared					                            # partition or queue name
#SBATCH --output=%j.out						                            # file to collect standard output
#SBATCH --error=%j.err						                            # file to collect standard output
#SBATCH --job-name=EmptyDrops					                        # name job for easier spotting, controlling
#SBATCH --mem=30GB						                                # each job from the array will get its own private 30G to work with)

# You may not place any commands before the last SBATCH directive

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
module load R

# list modules
module list

echo "== This is the script =="
Rscript 02_empty_droplets_stats_loop.R $id
echo "== End of Job =="

## This script was made using JHPCE 3.0 SLURM Cluster
## Available from http://xxxx
## CSC. Aug 14th, 2023

