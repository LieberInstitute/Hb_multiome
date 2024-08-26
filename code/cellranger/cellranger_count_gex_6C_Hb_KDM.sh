#!/bin/bash
#SBATCH --partition=katun		              # partition or queue name
#SBATCH --output=S4-CR_count-%j.out	        # file to collect standard output
#SBATCH --error=S4-CR_count-%j.err	        # file to collect standard output
#SBATCH --job-name=S4_CR_count	            # name job for easier spotting, controlling
#SBATCH --cpus-per-task=4		            # number of cores
#SBATCH --mem=20GB			                # memory per __node__
#SBATCH --mem=80GB			                # memory per __node__

# You may not place any commands before the last SBATCH directive

# NOTE: This script runs IN THE SAME DIRECTORY in which you ran sbatch
#       So include a cd command to ensure that you run it in the expected
#	directory (where your data files are located) OR use absolute paths
#	when specifying your input and output files

# For troubleshooting, it is useful to include    date;hostname;pwd
# at the start of your script. This sample script just runs those cmds.

echo "**** Job starts ****"
SAMPLE=S4C_Hb_KDM
echo "Run cellranger-count for only GEX side from a multiome experiment"
echo "https://kb.10xgenomics.com/hc/en-us/articles/360059656912-Can-I-analyze-only-the-Gene-Expression-data-from-my-single-cell-multiome-experiment"
echo ""
echo "This jobs is setup to run with 4 cores (20G mem) to get 80G mem"
echo "Sample: ${SAMPLE} 2024 Human Habenula"
echo "Total dataset size: xG GEX"
echo "Habenula folder name: ${SAMPLE}"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## load CellRanger-ATAC
module load cellranger/7.2.0

## List current modules for reproducibility
module list

## Locate file
echo "Processing sample ${SAMPLE}"
mkdir -p /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellranger/${SAMPLE}/
date

## Run cellRanger count
cellranger count --id=${SAMPLE} \
    --reference=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
    --fastqs=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024/GEX/${SAMPLE}/ \
    --chemistry=ARC-v1 \
    --localcores=4 \
    --localmem=80

## Move output
echo "Moving data to new location"
#mkdir -p /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/processed-data/cellranger_run_fast_version/
mv ${SAMPLE} /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/cellrangerGEX/
echo "Data moved to new location"

echo "**** Job ends ****"
date

## This script was made for slurm
## For running:
##      $ sbatch cellranger_count_gex_4C_Hb_KDM.sh
## CSC Aug, 2024
