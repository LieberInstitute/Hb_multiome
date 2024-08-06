#!/bin/bash
#SBATCH --partition=shared		            # partition or queue name
# SBATCH --time=1:00:00			            # time limit (D-HH:MM:SS)
#SBATCH --output=S2-cellranger-%j.out		# file to collect standard output
#SBATCH --error=S2-cellranger-%j.err		# file to collect standard output
#SBATCH --job-name=S2_cellrangerarc	        # name job for easier spotting, controlling
#SBATCH --cpus-per-task=8		            # number of cores
# SBATCH --ntasks=1			                # number of tasks running in parallel
# SBATCH --nodes=1			                # number of nodes
#SBATCH --mem=128GB			                # memory per __node__

# You may not place any commands before the last SBATCH directive

# NOTE: This script runs IN THE SAME DIRECTORY in which you ran sbatch
#       So include a cd command to ensure that you run it in the expected
#	directory (where your data files are located) OR use absolute paths
#	when specifying your input and output files

# For troubleshooting, it is useful to include    date;hostname;pwd
# at the start of your script. This sample script just runs those cmds.

echo "**** Job starts ****"
echo "This jobs is setup to run with 8 cores (16G mem) to get 128G mem"
echo "Sample: S2_Hb_KDM Human Habenula"
echo "Total dataset size: 24G ATAC / 21G GEX"
echo "HPC folder name: S2_Hb_KDM"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
#echo "Task id: ${$SLURM_ARRAY_TASK_ID}"

## load CellRanger
module load cellranger_arc/2.0.2

## List current modules for reproducibility
module list

## Locate file
SAMPLE=S2_Hb_KDM
echo "Processing sample ${SAMPLE}"
mkdir -p /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellranger/${SAMPLE}/
#mkdir -p /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellranger/${SAMPLE}/logs/
echo "Log folder created"
date

## Run CellRanger
cellranger-arc count --id=${SAMPLE} \
    --reference=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
    --libraries=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellranger/libraries_${SAMPLE}.csv \
    --localcores=8 \
    --localmem=128

## Move output
echo "Moving data to new location"
#mkdir -p /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/processed-data/cellranger_run_fast_version/
mv ${SAMPLE} /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/cellrangerARC/
echo "Data moved to new location"

echo "**** Job ends ****"
date

## This script was made for slurm version xxx
## CSC Jan 15, 2024
