#!/bin/bash
#SBATCH --partition=katun		              # partition or queue name
#SBATCH --job-name=CellR_ARC	            # name job for easier spotting, controlling
#SBATCH --cpus-per-task=8		              # number of cores
#SBATCH --mem=128GB			                  # memory per __node__
#SBATCH --array=1-2                       #change the number 2 to the number of entries in array_targets.txt
SAMPLE=$( sed -n ${SLURM_ARRAY_TASK_ID}p cellranger_arc_samples2024_package2_S7_S8.txt)   
#SBATCH --output=/logs/CellR_ARC-%j.out	        # file to collect standard output
#SBATCH --error=/logs/CellR_ARC-%j.err	      # file to collect standard output

# NOTE: This script runs IN THE SAME DIRECTORY in which you ran sbatch
#       So include a cd command to ensure that you run it in the expected
#	directory (where your data files are located) OR use absolute paths
#	when specifying your input and output files

# For troubleshooting, it is useful to include    date;hostname;pwd
# at the start of your script. This sample script just runs those cmds.

echo "**** Job starts ****"
echo "This script runs cellranger-arc v2"
echo "Processing cellranger-ARC experiments"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job name (script): ${SLURM_JOB_NAME}"
echo "Hostname (computer node): ${HOSTNAME}"
echo ""
echo "Job id: ${SLURM_JOBID}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "Array/sample: $SAMPLE"

## load CellRanger-ARC
module load cellranger_arc/2.0.2
module list


echo "This jobs is setup to run with 8 cores (16G mem) to get 128G mem"
echo "Sample: ${SAMPLE} 2024 package 2 / Human Habenula"
date
## Locate file
echo "Processing sample ${SAMPLE}"
mkdir -p /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellranger/${SAMPLE}/
echo "Habenula directory: ${SAMPLE}"
#echo "Log folder created"

### Run CellRanger-arc
cellranger-arc count --id=${SAMPLE} \
    --reference=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
    --libraries=/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellranger/multiome_library_${SAMPLE}.csv \
    --localcores=8 \
    --localmem=128

### Move output
echo "Moving data to new location"
mv ${SAMPLE} /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/cellrangerARC/
echo "Data moved to new location"

echo "**** Job ends ****"
date

## This script was made for slurm 22.05.9
## CSC Sep., 2024
