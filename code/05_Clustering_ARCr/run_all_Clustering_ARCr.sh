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
module load conda_R/4.3.x
## List current modules for reproducibility
module list

MAINDIR="/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome"
CODEDIR="${MAINDIR}/code"
PROCESSEDIR="${MAINDIR}/processed-data"
PLOTDIR="${MAINDIR}/plots"

## Update code style
# cd ${CODEDIR}
# Rscript update_style.R

## change dir
SUBDIR="05_Clustering_ARCr"
cd ${CODEDIR}/${SUBDIR}

echo "===== Run clustering with different methods ########"
echo "Current dir:"
pwd
## rm previous log files and output files
rm -f ${CODEDIR}/${SUBDIR}/logs/01_clustering_std_method_*.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/01_clustering_std_method/*.rds
rm -f ${PROCESSEDIR}/${SUBDIR}/01_clustering_std_method/cvs_files/*.csv
rm -f ${PLOTDIR}/${SUBDIR}/01_clustering_std_method/*.png
echo "Previous logs, cvs and png files deleted!"

## Run independent job / scritpt 'v2' compute only 24 (a subset) WNN configurations
sbatch 01_clustering_std_method_v2.sh   
# id1=$(sbatch --parsable 01_clustering_std_method.sh)      # perform 48 WNN configurations
# id1=$(sbatch --parsable 01_clustering_std_method_v2.sh)   # this 'v2' compute only 24 (a subset) WNN configurations
# echo $id1

echo "Clustering done!!"


echo "===== Running cell-types identification ..................."

# ===== Extended version to identy cell-types from different WNN methods (Louvain, LouvainM, SML and Leiden) 
# ===== Requiered for EDA (Louvain, LouvainM, SML and Leiden) 
#
# ## rm previous log files and output files
# rm ${CODEDIR}/${SUBDIR}/logs/02_v2_Hb_celltypes_ARCr*.txt
# rm ${PROCESSEDIR}/${SUBDIR}/02_Hb_celltypes_from_seurat_reanalyze_v2/*.csv
# rm ${PROCESSEDIR}/${SUBDIR}/02_Hb_celltypes_from_seurat_reanalyze_v2/cvs_files_markers/*.csv
# echo "Previous logs and output files deleted!"
# 
# # Run dependency job
# sbatch 02_Hb_celltypes_from_seurat_reanalyze_v2.sh
# echo "Cell type identification v2 done!!"


# ===== Version to identy cell-types only for WNN methods with Leiden

## rm previous log files and output files
rm -f ${CODEDIR}/${SUBDIR}/logs/02_Hb_celltypes_ARCr_v3*.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/02_Hb_celltypes_from_seurat_reanalyze_v3/cvs_files_markers/*.csv
echo "Previous logs and output files deleted!"

# Run dependency job
sbatch 02_Hb_celltypes_from_seurat_reanalyze_v3.sh
echo "Cell type identification v3 done!!"


echo "===== Running cell-types percentages ..................."

## rm previous log files and output files
## rm previous log files and output files
rm -f ${CODEDIR}/${SUBDIR}/logs/02_Hb_celltypes_ARCr_v3*.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/02_Hb_celltypes_from_seurat_reanalyze_v3/cvs_files_markers/DETAIL*v3.csv
rm -f ${PROCESSEDIR}/${SUBDIR}/02_Hb_celltypes_from_seurat_reanalyze_v3/FULL_SUMMARY*v3.csv
echo "Previous logs and output files deleted!"

# Run dependency job
sbatch 03_cell_types_perc_WNN_leiden_harmony_v3.sh
echo "Cell type percentages for v3 done!!"

# ===== Extended version to identy cell-types from different Cell Ranger Pipelines. It includes:
# ===== Requiered for EDA
# cellranger_pipe=="CR_crossBarcodes"
# cellranger_pipe=="CR_complementBarcodes"
# cellranger_pipe=="CR_arc_reanalyze"
# cellranger_pipe=="CR_arc_reanalyze_outliers"
# 03_cell_types_percentages.R

# ===== Requiered for EDA: without annotated cell-types
# 03_jaccard.sh
# ===== Requiered for EDA: without annotated cell-types /  called different outputs: Louvain, LouvainM, SML and Leiden
# 04_wnn_gene_expression.R


echo "===== Renaming idents ..................."

SUBsubDIR="05_rename_idents"

## rm previous log files and output files
rm -f ${CODEDIR}/${SUBDIR}/logs/05_rename_idents_*.txt
rm -f ${PROCESSEDIR}/${SUBDIR}/${SUBsubDIR}/*.rds
rm -f ${PLOTDIR}/${SUBDIR}/${SUBsubDIR}/*.pdf
echo "Previous logs and output files deleted!"

# Run dependency job
id3=$(sbatch --parsable 05_rename_idents.sh)
# Alternative command if you need to re-run the clustering
# id3=$(sbatch --parsable --dependency=afterok:$id2 05_rename_idents.sh)
echo $id3
echo "Renaming Idents Done!!"


echo "===== Processing Silhouette and Jaccard ..................."

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


echo "===== Processing visualization plots on renamed WNN idents curated  ..................."
#Generate WNN VPlots, FeaturePlots, Other Jaccard, etc on renamed Idents curated (Visium-Project)

SUBsubDIR="08_wnn_gene_expression_plts_renamed_idents"

## rm previous log files and output files
# rm ${CODEDIR}/${SUBDIR}/logs/06_jaccard_annotated_*.txt
rm ${PROCESSEDIR}/${SUBDIR}/${SUBsubDIR}/Silhouette_*.cvs
rm ${PLOTDIR}/${SUBDIR}/${SUBsubDIR}/*.pdf
rm ${PLOTDIR}/${SUBDIR}/${SUBsubDIR}/*.png
echo "Previous logs and output files deleted!"

# Run dependency job
Rscript 08_wnn_gene_expression_plts_renamed_idents.R
# id4=$(sbatch --parsable --dependency=afterok:$id3 08_wnn_gene_expression_plts_renamed_idents.sh
echo "08_wnn_gene_expression_plts_renamed_idents done!!"


echo "===== Processing 10_format_DEG_csv_files.R - rename columns in top50-DEG CSV files for better clarity and consistency"

SUBsubDIR="02_Hb_celltypes_from_seurat_reanalyze_v3/cvs_files_formatted"
echo ${PROCESSEDIR}/${SUBDIR}/${SUBsubDIR}
## rm previous log files and output files
rm -f ${PROCESSEDIR}/${SUBDIR}/${SUBsubDIR}/WNN*top50.csv
echo "Previous logs and output files deleted!"

sbatch 10_format_DEG_csv_files.sh
echo "10_format_DEG_csv_files.R Done!!"


echo "**** Job ends ****"
date


} > $log_path 2>&1



## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/