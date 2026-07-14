#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=2G
#SBATCH --job-name=03_link_inputs
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/16_shiny_app/logs/03_link_inputs.txt
#SBATCH -e ../../processed-data/16_shiny_app/logs/03_link_inputs.txt

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
dest_dir=${repo_dir}/code/16_shiny_app

#   Using relative paths for sym links so things work with git
cd ${dest_dir}

f_name=merged_metacell_seur.qs2
rm -f ${f_name}
ln -s ../../processed-data/16_shiny_app/01_prep_objects/${f_name} ${f_name}

f_name=trios.csv
rm -f ${f_name}
ln -s ../../processed-data/16_shiny_app/01_prep_objects/${f_name} ${f_name}

f_name=atlas_seur_minimal.qs2
rm -f ${f_name}
ln -s ../../processed-data/16_shiny_app/02_prep_cell_level/${f_name} ${f_name}

echo "**** Job ends ****"
date
