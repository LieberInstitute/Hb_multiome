#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=04_search_peaks
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL

log_path=logs/04_search_peaks.txt

{
set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module (used for Seurat)
module load conda_R/4.3.x

# In local env use your user library where Seurat was installed & source & activate your macs2 custom env
# As I am working with Modules in shared env: to avoid conflicts with shared and custom conda envs, 
#    I am calling macs2 directly from CallPeaks( function() into the R script

# if local =====================================================================
# export R_LIBS_USER="/users/csoto/R/4.3.x"
#eval "$(conda shell.bash hook)"
## Activate Conda for MACS2
# source /home/csoto/miniconda3/etc/profile.d/conda.sh
# conda activate macs2_conda3_env

# Confirm MACS2
# echo "Conda env: $CONDA_DEFAULT_ENV"
# which macs2
# macs2 --version
# if local =====================================================================

## Reproducibility info
module list

which /users/csoto/.conda/envs/macs2_conda3_env/bin/macs2
ls /users/csoto/.conda/envs/macs2_conda3_env/bin/macs2

## Your main script
Rscript 04_search_peaks.R

# Capture return code and exit safely
ret=$?
echo "**** Job ends ****"
date
echo "Exit code: $ret"
exit $ret

} > $log_path 2>&1
