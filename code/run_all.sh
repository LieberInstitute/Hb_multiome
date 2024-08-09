#!/bin/bash
#SBATCH --partition=katun	
#SBATCH --job-name=Hb_allQCs
# SBATCH -c 1
# SBATCH -t 1-00:00:00
#SBATCH --mem=30GB						                                    # each job from the array will get its own private 30G to work with)
#SBATCH --array=1-5                                           	  # change the array range acordingly with the rows listed in the array_targets.txt
id=$( sed -n ${SLURM_ARRAY_TASK_ID}p ../array_targets_names.txt)				                                  
#SBATCH --output=logs/Hb_allQC%j.out						                              
#SBATCH --error=logs/Hb_allQC_%j.err		
#SBATCH --mail-type=ALL

echo "**** Job starts ****"
echo "General script to QCed CellRanger-ARC datasets"
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
PLOTDIR="${MAINDIR}/plots"


## Update code style
# cd ${CODEDIR}
# Rscript update_style.R


## Analysis workflow for `Multiome 10x CellRanger scRNAseq + scATACseq data`

## First we calculate estimated number of cells by EmptyDroplets for comparison purposes 
## Delete the logs/old-results, and re-submit EmptyDrops
cd ${CODEDIR}/00_RNAbackground
# mkdir -p logs ## Create the logs directory if it doesn't exist
rm logs/*.err
rm logs/*.out
rm ${PROCESSEDIR}/00_RNAbackground/*.csv
rm ${PLOTDIR}/00_RNAbackground/*.png
## This is an independent jobs-array 
sbatch 01_empty_drops_in_multiomedata.sh


# ## Add future steps here
# 
# ## Read in some raw-data into R
# cd ${CODEDIR}/01_read_data_to_r
# mkdir -p logs ## Create the logs directory if it doesn't exist
# # Delete output data before re-generating them:
# rm -f ${PROCESSEDIR}/01_read_data_to_r/penguins_data.csv
# id1=$(sbatch --parsable 01_read_data_to_r.sh)
# 
# ## Explore the data
# cd ${CODEDIR}/02_explore_data
# mkdir -p logs ## Create the logs directory if it doesn't exist
# id2=$(sbatch --parsable --dependency=afterok:$id1 01_ggpairs.sh)
# sbatch --dependency=afterok:$id2 02_boxplots.sh

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/

