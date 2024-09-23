#!/bin/bash
#SBATCH --partition=katun		              
#SBATCH --job-name=cellranger_arc_reanalyze_S3_Hb
#SBATCH --cpus-per-task=8		              # number of cores
#SBATCH --mem=128GB			                  # memory per __node__
#SBATCH --output=/logs/cellranger_arc_reanalyze_S3_Hb_%j.txt
#SBATCH --error=/logs/cellranger_arc_reanalyze_S3_Hb_%j.txt

# This jobs is setup to run with 8 cores (16G mem) to get 128G mem

echo "**** Job starts ****"
SAMPLE=S3_Hb_KDM
BC=3C_Hb_KDM
echo "Sample: ${SAMPLE} 2024 Human Habenula"
echo "Habenula folder name: ${SAMPLE}"
date


echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## load CellRanger-ARC
module load cellranger_arc/2.0.2
module list

## Locate file
echo "Processing sample ${SAMPLE}"
# mkdir -p /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellrangerARC_reanalyze/${SAMPLE}/
# echo "Log folder created"
date

## Run CellRanger
cellranger-arc reanalyze --id=${SAMPLE}_reanalysis \
    --barcodes=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/100_cell_match_cellranger_gex_atac/reanalize_files/${BC}_matching_ARC_barcodes_1.csv \
    --matrix=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/cellrangerARC/${SAMPLE}/outs/raw_feature_bc_matrix.h5 \
    --atac-fragments=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/cellrangerARC/${SAMPLE}/outs/atac_fragments.tsv.gz \
    --reference=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
    --localcores=8 \
    --localmem=128

## Move output
echo "Moving data to new location"
mv ${SAMPLE}_reanalysis /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/cellrangerARC_reanalyze/
echo "Data moved to new location"

echo "**** Job ends ****"
date

## This script was made for slurm 22.05.9
## CSC Sep., 2024
