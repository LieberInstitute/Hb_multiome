#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=10G
#SBATCH --job-name=05_gene_set_analysis
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-15%10

#   Run just the gene-set-analysis step from MAGMA for each GWAS and cell-type
#   resolution

## Define loops and appropriately subset each variable for the array task ID
all_gwas=(MDD MDD2019 panic SCZ SUD2020)
gwas=${all_gwas[$(( $SLURM_ARRAY_TASK_ID / 3 % 5 ))]}

all_cell_type_group=(broad semi_broad mid)
cell_type_group=${all_cell_type_group[$(( $SLURM_ARRAY_TASK_ID / 1 % 3 ))]}

## Explicitly pipe script output to a log
log_path=../../processed-data/10_MAGMA/logs/05_gene_set_analysis_${gwas}_${cell_type_group}_${SLURM_ARRAY_TASK_ID}.txt

{
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
out_prefix=${repo_dir}/processed-data/10_MAGMA/$gwas/$cell_type_group
gene_set_path=${repo_dir}/processed-data/10_MAGMA/gene_sets/${cell_type_group}.tsv
gene_results_path=${repo_dir}/processed-data/10_MAGMA/$gwas/${gwas}.genes.raw

module load magma/1.10

## List current modules for reproducibility
module list

#   Gene set analysis step
magma \
    --gene-results $gene_results_path \
    --set-annot $gene_set_path gene-col=gene_id set-col=set_id \
    --out $out_prefix

echo "**** Job ends ****"
date

} > $log_path 2>&1

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/
