#!/bin/bash
#SBATCH --partition=katun		              # partition or queue name
#SBATCH --job-name=CellR_atac	            # name job for easier spotting, controlling
#SBATCH --cpus-per-task=8		              # number of cores
#SBATCH --mem=128GB			                  # memory per __node__
#SBATCH --array=1-6                       #change the number of entries
SAMPLE=$( sed -n ${SLURM_ARRAY_TASK_ID}p cellranger_atac_S7-S12.txt)   
#SBATCH --output=/logs/CellR_atac-%j.out	    
#SBATCH --error=/logs/CellR_atac-%j.err	      

# You may not place any commands before the last SBATCH directive

# NOTE: This script runs IN THE SAME DIRECTORY in which you ran sbatch
#       So include a cd command to ensure that you run it in the expected
#	directory (where your data files are located) OR use absolute paths
#	when specifying your input and output files

# For troubleshooting, it is useful to include    date;hostname;pwd
# at the start of your script. This sample script just runs those cmds.

echo "**** Job starts ****"
echo "This script runs cellranger-atac count in S7 to S12 (data-package2)"
echo "Processing cellranger-atac for multiome atac data"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job name (script): ${SLURM_JOB_NAME}"
echo "Hostname (computer node): ${HOSTNAME}"
echo ""
echo "Job id: ${SLURM_JOBID}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "Array/sample: $SAMPLE"

## load CellRanger-count
module load cellranger-atac/2.1.0
module list

## Locate file
echo "Processing sample ${SAMPLE}"
mkdir -p /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellranger/${SAMPLE}/
date

echo "Running cellranger count for sample ${SAMPLE}"
## Run CellRanger-ATAC
cellranger-atac count --id=${SAMPLE} \
    --reference=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
    --fastqs=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024_data_package2/ATAC/${SAMPLE} \
    --chemistry=ARC-v1 \
    --localcores=4 \
    --localmem=80

## Move output
echo "Moving data to new location"
mv ${SAMPLE} /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/cellrangerATAC/
echo "Data moved to new location: `~/processed-data/cellrangerATAC/`"

echo "**** Job ends ****"
date

## This script was made for slurm
## For running:
##      $ sbatch cellranger_count_atac_S7-S12.sh
## CSC Sep, 2024
