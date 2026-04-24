#!/bin/bash
#SBATCH -p katun
#SBATCH --mem=20G
#SBATCH --job-name=05_mean_ratio_MAGMA
#SBATCH -c 1
#SBATCH -t 1-0:00:00
#SBATCH -o ../../../processed-data/10_MAGMA/RNA/logs/05_mean_ratio_MAGMA_%a.txt
#SBATCH -e ../../../processed-data/10_MAGMA/RNA/logs/05_mean_ratio_MAGMA_%a.txt
#SBATCH --array=1-16%16
#SBATCH --reservation=neagles-2wk

#   Run all 3 steps in the MAGMA pipeline for every GWAS. Critically,
#   all relevant inputs/ reference files use hg19 and European ancestry.

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
hb_gwas_dir=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/10_MAGMA/habenula_pilot_gwas
multiome_dir=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome
out_dir=${repo_dir}/processed-data/10_MAGMA/RNA/$gwas
gene_loc=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/no_secondary/MAGMA/hg19_gene_loc.tsv
bfile=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/10_MAGMA/habenula_pilot_gwas/g1000_eur/g1000_eur

case ${gwas} in
    MDD)
        snp_loc=${hb_gwas_dir}/mdd2019edinburgh/PGC_UKB_depression_genome-wide.snploc
        pval_file=${multiome_dir}/processed-data/10_MAGMA/MDD2019/p_values.tsv
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
        #   Other GWASes follow a simpler pattern
        snp_loc=${multiome_dir}/processed-data/10_MAGMA/${gwas}/SNPs.tsv
        pval_file=${multiome_dir}/processed-data/10_MAGMA/${gwas}/p_values.tsv
        ;;
esac

echo "Processing GWAS ${gwas}"
echo "Using SNP location file: ${snp_loc}"
echo "Using p-value file: ${pval_file}"

mkdir -p $out_dir
mkdir -p $out_dir/mean_ratio

#   Annotation step
magma \
    --annotate \
    --snp-loc $snp_loc \
    --gene-loc $gene_loc \
    --out $out_dir/mean_ratio_$gwas

#   Gene analysis step
magma \
    --bfile $bfile \
    --pval $pval_file ncol=N \
    --gene-annot $out_dir/$gwas.genes.annot \
    --out $out_dir/mean_ratio_$gwas

#   Gene set analysis step

for i in broad mid fine; do
    echo "Running gene set analysis set for $i resolution"

    gene_set_path=${repo_dir}/processed-data/10_MAGMA/RNA/gene_sets/$i.tsv

    magma \
        --gene-results $out_dir/${gwas}.genes.raw \
        --set-annot $gene_set_path gene-col=gene_id set-col=set_id \
        --out $out_dir/mean_ratio/$i
done

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.4
## available from http://research.libd.org/slurmjobs/
