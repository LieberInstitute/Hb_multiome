#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=25G
#SBATCH --job-name=04_OUD_MAGMA
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/10_MAGMA/logs/04_OUD_MAGMA.txt
#SBATCH -e ../../processed-data/10_MAGMA/logs/04_OUD_MAGMA.txt

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

## List current modules for reproducibility
module list

repo_dir=$(git rev-parse --show-toplevel)
out_dir=${repo_dir}/processed-data/10_MAGMA/OUD

snp_loc=$repo_dir/processed-data/10_MAGMA/OUD/SNPs.tsv
pval_file=$repo_dir/processed-data/10_MAGMA/OUD/p_values.tsv
gene_loc=/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/geneloc/GRCh38_Ensembl-93_GENES_chr-x-y-mt.gene.loc
bfile=/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur

mkdir -p $out_dir

#   Annotation step
magma \
    --annotate \
    --snp-loc $snp_loc \
    --gene-loc $gene_loc \
    --out $out_dir/OUD

#   Gene analysis step
magma \
    --bfile $bfile \
    --pval $pval_file ncol=N \
    --gene-annot $out_dir/OUD.genes.annot \
    --out $out_dir/OUD

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/
