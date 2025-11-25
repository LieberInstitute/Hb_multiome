#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=25G
#SBATCH --job-name=04_MAGMA_first_two
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../processed-data/10_MAGMA/logs/04_MAGMA_first_two_%a.txt
#SBATCH -e ../../processed-data/10_MAGMA/logs/04_MAGMA_first_two_%a.txt
#SBATCH --array=1-5%5
#SBATCH --exclude=compute-175

#   Run the first two steps in the MAGMA pipeline for every GWAS. Critically,
#   all relevant inputs/ reference files use hg19 and European ancestry. This
#   needs to be rerun (instead of grabbing what was done for habenula pilot)
#   because of two issues: a lack of chr22 in the geneloc input and use of an
#   hg38 input file for the MDD GWAS

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

all_gwas=(MDD MDD2019 panic SCZ SUD2020)
gwas=${all_gwas[$(($SLURM_ARRAY_TASK_ID - 1))]}

repo_dir=$(git rev-parse --show-toplevel)
hb_gwas_dir=/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/13_MAGMA/GWAS
out_dir=${repo_dir}/processed-data/10_MAGMA/$gwas
gene_loc=${repo_dir}/processed-data/10_MAGMA/hg19_gene_loc.tsv
bfile=/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur

case ${gwas} in
    MDD)
        snp_loc=${repo_dir}/processed-data/10_MAGMA/MDD/SNPs.tsv
        pval_file=${repo_dir}/processed-data/10_MAGMA/MDD/p_values.tsv
        ;;
    MDD2019)
        snp_loc=${hb_gwas_dir}/mdd2019edinburgh/PGC_UKB_depression_genome-wide.snploc
        pval_file=${repo_dir}/processed-data/10_MAGMA/MDD2019/p_values.tsv
        ;;
    panic)
        snp_loc=${hb_gwas_dir}/panic2019/pgc-panic2019.snploc
        pval_file=${hb_gwas_dir}/panic2019/pgc-panic2019.pval
        ;;
    SCZ)
        snp_loc=${hb_gwas_dir}/scz2022/PGC3_SCZ_wave3.european.autosome.public.v3.snploc
        pval_file=${hb_gwas_dir}/scz2022/PGC3_SCZ_wave3.european.autosome.public.v3.pval
        ;;
    SUD2020)
        snp_loc=${hb_gwas_dir}/sud2020op/opi.DEPvEXP_EUR.noAF.snploc
        pval_file=${hb_gwas_dir}/sud2020op/opi.DEPvEXP_EUR.noAF.pval
        ;;
    *)
        #   Unexpected values
        echo "Unknown GWAS: ${gwas}"
        exit 1
        ;;
esac

echo "Processing GWAS ${gwas}"
echo "Using SNP location file: ${snp_loc}"
echo "Using p-value file: ${pval_file}"

mkdir -p $out_dir

#   Annotation step
magma \
    --annotate \
    --snp-loc $snp_loc \
    --gene-loc $gene_loc \
    --out $out_dir/$gwas

#   Gene analysis step
magma \
    --bfile $bfile \
    --pval $pval_file ncol=N \
    --gene-annot $out_dir/$gwas.genes.annot \
    --out $out_dir/$gwas

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/
