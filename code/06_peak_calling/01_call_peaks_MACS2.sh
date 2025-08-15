#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=80G
#SBATCH --job-name=01_call_peaks_MACS2
#SBATCH -c 2
#SBATCH -t 3-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL
#SBATCH --array=0-2%3   # 3 wnn resolution levels

set -eo pipefail

wnn_resolution=(Fine Broad Mid)

i=${SLURM_ARRAY_TASK_ID}
m=${#wnn_resolution[@]}

# guard
if (( i < 0 || i >= m )); then
  echo "Invalid SLURM_ARRAY_TASK_ID=$i (must be 0..$((m-1)))"
  exit 1
fi

res="${wnn_resolution[$i]}"

log_path=logs/01_call_peaks_MACS2_${res}_task_${i}.txt

{

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
#ls /users/csoto/.conda/envs/macs2_conda3_env/bin/macs2

## Your main script
Rscript 01_call_peaks_MACS2.R --wnn_resolution "${res}"

# Capture return code and exit safely
ret=$?
echo "**** Job ends ****"
date
echo "Exit code: $ret"
exit $ret

} > $log_path 2>&1
