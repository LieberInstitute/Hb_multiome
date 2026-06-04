#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=20G
#SBATCH --job-name=01_run_MAGMA
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_MAGMA/trios/logs/01_run_MAGMA_%a.txt
#SBATCH -e ../../../processed-data/10_MAGMA/trios/logs/01_run_MAGMA_%a.txt
#SBATCH --array=1-16%16

#   Run just the gene-set analysis step from MAGMA (other steps computed already
#   from other analyses)

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load magma/1.10
module list

all_gwas=(MDD panic SCZ SUD2020 AUD CUD ext_cannabis lifetime_cannabis OUD SUD2 compulsive internalizing neurodev p_factor SCZ_BPD SUD3)
gwas=${all_gwas[$(($SLURM_ARRAY_TASK_ID - 1))]}

repo_dir=$(git rev-parse --show-toplevel)
in_path=${repo_dir}/processed-data/10_MAGMA/RNA/$gwas/${gwas}.genes.raw
out_dir=${repo_dir}/processed-data/10_MAGMA/trios/$gwas

echo "Processing GWAS ${gwas}"

mkdir -p $out_dir

for resolution in broad fine; do
    echo "Running gene set analysis set for $resolution resolution"

    gene_set_path=${repo_dir}/processed-data/13_tripod_trios/11_GO/gene_sets_${resolution}.tsv

    magma \
        --gene-results $in_path \
        --set-annot $gene_set_path gene-col=gene_id set-col=set_id \
        --out $out_dir/$resolution
done

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/
