#!/bin/bash
#SBATCH --partition=shared		          # partition or queue name
#SBATCH --output=cellranger-A%j.out	    # file to collect standard output
#SBATCH --error=cellranger-A%j.err	    # file to collect standard output
#SBATCH --job-name=cellrangerarc-A	    # name job for easier spotting, controlling
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
SAMPLE=3A_Hb_KDM
echo "This jobs is setup to run with 4 cores (20G mem) to get 80G mem"
echo "Sample: ${SAMPLE} 2024 Human Habenula"
echo "Total dataset size: 34G ATAC"
echo "Habenula folder name: ${SAMPLE}"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## load CellRanger-ATAC
module load cellranger-atac/2.1.0

## List current modules for reproducibility
module list

## Locate file
echo "Processing sample ${SAMPLE}"
mkdir -p /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellranger/${SAMPLE}/
date

## Run CellRanger-ATAC
cellranger-atac count --id=${SAMPLE} \
    --reference=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
    --fastqs=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024/ATAC/${SAMPLE}/ \
    --chemistry=ARC-v1 \
    --localcores=4 \
    --localmem=80

## Move output
echo "Moving data to new location"
#mkdir -p /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/processed-data/cellranger_run_fast_version/
mv ${SAMPLE} /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/cellrangerATAC/
echo "Data moved to new location"

echo "**** Job ends ****"
date

## This script was made for slurm
## For running:
##      $ sbatch <file_name.sh>
## CSC Aug 05, 2024
