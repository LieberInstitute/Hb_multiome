#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=25G
#SBATCH --job-name=01_first_test_template
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-21%20

#   Using pre-computed gene-level statistics from the habenula pilot project,
#   just run the gene-set analysis step from MAGMA

## Define loops and appropriately subset each variable for the array task ID
all_gwas=(MDD MDD2019 OUD panic SCZ SUD SUD2020)
gwas=${all_gwas[$(( $SLURM_ARRAY_TASK_ID / 3 % 7 ))]}

all_cell_type_group=(broad semi_broad mid)
cell_type_group=${all_cell_type_group[$(( $SLURM_ARRAY_TASK_ID / 1 % 3 ))]}

## Explicitly pipe script output to a log
log_path=../../processed-data/10_MAGMA/logs/01_first_test_template_${gwas}_${cell_type_group}_${SLURM_ARRAY_TASK_ID}.txt

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

hb_gwas_dir=/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS
out_prefix=${repo_dir}/processed-data/10_MAGMA/$gwas/$cell_type_group
gene_set_path=${repo_dir}/processed-data/10_MAGMA/gene_sets/${cell_type_group}.tsv

mkdir -p $(dirname $out_prefix)

#   Processed gene-level statistics from MAGMA in the habenula pilot project
case ${gwas} in
    MDD)
        gene_results_path=${hb_gwas_dir}/MDD/MDD.phs001672.pha005122.genes.raw
        ;;
    MDD2019)
        gene_results_path=${hb_gwas_dir}/mdd2019edinburgh/PGC_UKB_depression_genome-wide.genes.raw
        ;;
    OUD)
        #   TODO: don't have path yet
        ;;
    panic)
        gene_results_path=${hb_gwas_dir}/panic2019/pgc-panic2019.genes.raw
        ;;
    SCZ)
        gene_results_path=${hb_gwas_dir}/scz2022/PGC3_SCZ_wave3.european.autosome.public.v3.ensembl.genes.raw
        ;;
    SUD)
        #   Might remove this
        ;;
    SUD2020)
        gene_results_path=${hb_gwas_dir}/sud2020op/opi.DEPvEXP_EUR.noAF.genes.raw
        ;;
    *)
        #   Unexpected values
        echo "Unknown GWAS: ${gwas}"
        exit 1
        ;;
esac

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

