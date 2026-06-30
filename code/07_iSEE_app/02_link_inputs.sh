#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=2G
#SBATCH --job-name=02_link_inputs
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/07_iSEE_app/logs/02_link_inputs.txt
#SBATCH -e ../../processed-data/07_iSEE_app/logs/02_link_inputs.txt

# Symlink Shiny app files so relative paths can be used from the destination
# directories for the apps

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

repo_dir=$(git rev-parse --show-toplevel)
rna_dest_dir=${repo_dir}/code/07_iSEE_app/RNA
atac_dest_dir=${repo_dir}/code/07_iSEE_app/ATAC
rna_f_name=sce_RNA_iSEE.qs2
atac_f_name=sce_ATAC_iSEE.qs2

#   Using relative paths for sym links so things work with git
cd ${rna_dest_dir}
rm -f ${rna_f_name}
ln -s ../../../processed-data/07_iSEE_app/01_prep_sce/${rna_f_name} ${rna_f_name}

cd ${atac_dest_dir}
rm -f ${atac_f_name}
ln -s ../../../processed-data/07_iSEE_app/01_prep_sce/${atac_f_name} ${atac_f_name}

echo "**** Job ends ****"
date
