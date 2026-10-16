#!/bin/bash
#SBATCH -p transfer
#SBATCH --mem=5G
#SBATCH --job-name=04_upload_raw
#SBATCH -c 1
#SBATCH -t 20-00:00:00
#SBATCH -o ../../processed-data/19_data_uploads/logs/04_upload_raw.txt
#SBATCH -e ../../processed-data/19_data_uploads/logs/04_upload_raw.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## List current modules for reproducibility
module list

repo_dir=$(git rev-parse --show-toplevel)

#   For authentication
source /users/neagles/private_vars.txt

ascp \
    -l 1000M \
    -k 2 \
    -QT \
    $repo_dir/processed-data/19_data_uploads/01_file_map/fastq_flat/*.fastq.gz \
    neagles7-restricted@aspera.nemoarchive.org:pt9rm3g/

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.4.0
## available from http://research.libd.org/slurmjobs/
