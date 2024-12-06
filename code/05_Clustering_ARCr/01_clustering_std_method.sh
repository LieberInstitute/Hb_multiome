#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=60G
#SBATCH --job-name=01_clustering_std_method
#SBATCH -c 2
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL
#SBATCH --array=1-48%20

## Define loops and appropriately subset each variable for the array task ID
all_clust_method=(1 2 3 4)
clust_method=${all_clust_method[$(( $SLURM_ARRAY_TASK_ID / 12 % 4 ))]}
# 1 = original Louvain algorithm
# 2 = Louvain algorithm with multilevel refinement
# 3 = SLM algorithm
# 4 = Leiden algorithm

all_clust_res=(0.8 1 1.5 2)
clust_res=${all_clust_res[$(( $SLURM_ARRAY_TASK_ID / 3 % 4 ))]}
# 0.8 = Default
# 1 = a bit larger clusters
# 1.5 =  more larger clusters
# 2 = more more larger clusters

all_knn=(20 30 40)
knn=${all_knn[$(( $SLURM_ARRAY_TASK_ID / 1 % 3 ))]}
# 20 = Default
# 30 = 30 NN
# 40 = 40 NN (Useful if the dataset is very noisy)

## Explicitly pipe script output to a log
log_path=logs/01_clustering_std_method_${clust_method}_${clust_res}_${knn}_${SLURM_ARRAY_TASK_ID}.txt

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

## Load the R module
module load conda_R/4.3.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_clustering_std_method.R --clust_method ${clust_method} --clust_res ${clust_res} --knn ${knn}

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/

