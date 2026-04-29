#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=150G
#SBATCH --job-name=03_tripod_trios
#SBATCH -c 8
#SBATCH -t 20-0:00:00
#SBATCH -o ../../processed-data/13_tripod_trios/logs/03_tripod_trios_%a.txt
#SBATCH -e ../../processed-data/13_tripod_trios/logs/03_tripod_trios_%a.txt
#SBATCH --array=18
#SBATCH --reservation=neagles-2wk

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5

## List current modules for reproducibility
module list

Rscript 03_tripod_trios.R

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/
#!/bin/bash