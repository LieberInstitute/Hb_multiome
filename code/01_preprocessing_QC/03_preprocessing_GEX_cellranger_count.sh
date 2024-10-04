#!/bin/bash -l
#SBATCH --partition=katun	
#SBATCH --job-name=03_preprocessing_GEX_cellranger_count
#SBATCH --mem=20GB						                            
# SBATCH --mail-user=bioinformatic2019@gmail.com
# SBATCH --mail-type=end,fail
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-12                                       #change the number x to the number of entries in array_targets.txt
id=$(sed -n ${SLURM_ARRAY_TASK_ID}p target_names_cellranger_count.txt)

## Explicitly pipe script output to a log
log_path=logs/03_preprocessing_GEX_cellranger_count_${SLURM_ARRAY_TASK_ID}.txt

echo "Removing previous outputs ..."
echo " "
## remove previous outputs if exists
MAINDIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome"
SUBDIR="01_preprocessing_QC/cellranger_count"
PROCESSEDIR="${MAINDIR}/processed-data/${SUBDIR}"
PLOTDIR="${MAINDIR}/plots/${SUBDIR}"


## Delete the logs/old-results, and re-submit seurat builder
# mkdir -p logs ## Create the logs directory if it doesn't exist
#rm logs/03_preprocessing_GEX_cellranger_count_*.txt
rm ${PROCESSEDIR}/*.csv
rm ${PROCESSEDIR}/*.rds
rm ${PLOTDIR}/*.png
echo " "


echo "**** Job starts ****"
echo "     Build Seurat objects and calculate basic statistics for GEX from multiome data run with cellranger-count"
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

Rscript 03_preprocessing_GEX_cellranger_count.R $id

echo "**** Job ends ****"
date


## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
## CSC. Oct 4th, 2024 
