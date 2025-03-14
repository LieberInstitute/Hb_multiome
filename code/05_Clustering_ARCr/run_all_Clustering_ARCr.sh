#!/bin/bash
#SBATCH --partition=katun	
#SBATCH --job-name=run_all_Clustering_ARCr
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=30GB
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL

## Explicitly pipe script output to a log
log_path=logs/05_rename_idents_${clustering_name}_${SLURM_ARRAY_TASK_ID}.txt


{
set -e

echo "**** Job starts ****"
echo "Run run_all_Clustering_ARCr"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job name (script): ${SLURM_JOB_NAME}"
echo "Hostname (computer node): ${HOSTNAME}"
echo ""
echo "Job id: ${SLURM_JOBID}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load modules
module conda_R/4.3.x
## List current modules for reproducibility
module list

MAINDIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome"
CODEDIR="${MAINDIR}/code"
PROCESSEDIR="${MAINDIR}/processed-data"
PLOTDIR="${MAINDIR}/plots"

## Update code style
# cd ${CODEDIR}
# Rscript update_style.R


######## Run clustering with different methods
######## All section runs on SUBDIR="05_Clustering_ARCr"

SUBDIR="05_Clustering_ARCr"

echo "Running clustering ..................."

# ## rm previous log files and output files
# rm ${CODEDIR}/${SUBDIR}/logs/01_clustering_std_method_*.txt
# rm ${PROCESSEDIR}/${SUBDIR}/01_clustering_std_method/*.rds
# rm ${PROCESSEDIR}/${SUBDIR}/01_clustering_std_method/cvs_files/*.csv
# rm ${PLOTDIR}/${SUBDIR}/01_clustering_std_method/*.png

# echo "Previous logs and output files deleted!"

## Run independent job

# id1=$(sbatch --parsable 01_clustering_std_method.sh)      # perform 48 WNN configurations
# id1=$(sbatch --parsable 01_clustering_std_method_v2.sh)   # this 'v2' compute only 24 (a subset) WNN configurations
# echo $id1

echo "Clustering done!!"


######## Run cell type identification

echo "Running cell-types identification ..................."

# ## rm previous log files and output files
# rm ${CODEDIR}/${SUBDIR}/logs/02_v2_Hb_celltypes_ARCr*.txt
# rm ${PROCESSEDIR}/${SUBDIR}/02_Hb_celltypes_from_seurat_reanalyze_v2/*.csv
# rm ${PROCESSEDIR}/${SUBDIR}/02_Hb_celltypes_from_seurat_reanalyze_v2/cvs_files_markers/*.csv
# 
# echo "Previous logs and output files deleted!"
# 
# # Run dependency job
# 
# id2=$(sbatch --parsable 02_Hb_celltypes_from_seurat_reanalyze_v2.sh)
# # Alternative command if you need to re-run the clustering
# # id2=$(sbatch --parsable --dependency=afterok:$id1 02_Hb_celltypes_from_seurat_reanalyze_v2.sh)
# 
# echo $id2

echo "Cell type identification v2 done!!"


######## Run rename idents

echo "Renaming idents ..................."

SUBsubDIR="05_rename_idents"

## rm previous log files and output files
rm ${CODEDIR}/${SUBDIR}/logs/05_rename_idents_*.txt
rm ${PROCESSEDIR}/${SUBDIR}/${SUBsubDIR}/*.rds
rm ${PLOTDIR}/${SUBDIR}/${SUBsubDIR}/*.pdf

echo "Previous logs and output files deleted!"

# Run dependency job

id3=$(sbatch --parsable 05_rename_idents.sh)
# Alternative command if you need to re-run the clustering
# id3=$(sbatch --parsable --dependency=afterok:$id2 05_rename_idents.sh)

echo $id3

echo "Renaming Idents Done!!"


######## Run Jaccard on renamed Idents for easy identification: 03_jaccard.sh

echo "Processing Silhouette and Jaccard ..................."

SUBsubDIR="06_jaccard_annotated"

## rm previous log files and output files
rm ${CODEDIR}/${SUBDIR}/logs/06_jaccard_annotated_*.txt
rm ${PROCESSEDIR}/${SUBDIR}/${SUBsubDIR}/Silhouette_*.cvs
rm ${PLOTDIR}/${SUBDIR}/${SUBsubDIR}/*.pdf
rm ${PLOTDIR}/${SUBDIR}/${SUBsubDIR}/*.png

echo "Previous logs and output files deleted!"

# Run dependency job

id4=$(sbatch --parsable --dependency=afterok:$id3 06_jaccard_annotated.sh)

echo $id4

echo "Silhouette and Jaccard done!!"


######## Run WNN GEX on renamed Idents for easy identification: 04_wnn_gene_expression.R


echo "**** Job ends ****"
date

} > $log_path 2>&1

## Script to renamed the columns on the DEG CVS files in a nice format

Rscript 10_format_DEG_csv_files.R

## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/