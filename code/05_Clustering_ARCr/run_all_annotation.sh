#!/bin/bash
#SBATCH --partition=katun	
#SBATCH --job-name=run_all_annotation
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=2GB
#SBATCH -o /dev/null
#SBATCH -e /dev/null
# SBATCH --mail-type=ALL

## Explicitly pipe script output to a log
log_path=logs/run_all_annotation.txt


{
set -e

echo "**** Job starts ****"
echo "Run run_all_annotation"
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

## change dir
SUBDIR="05_Clustering_ARCr"
cd ${CODEDIR}/${SUBDIR}

echo "===== Run clustering with different methods ########"
echo "Current dir:"
pwd
## rm previous log files and output files
rm -f ${CODEDIR}/${SUBDIR}/logs/16_rename_idents_post_mean_ratio_HClust.txt




## remove data objects
# rm -f ${PROCESSEDIR}/${SUBDIR}/01_clustering_std_method/*.rds
# rm -f ${PROCESSEDIR}/${SUBDIR}/01_clustering_std_method/cvs_files/*.csv
# rm -f ${PLOTDIR}/${SUBDIR}/01_clustering_std_method/*.png
echo "Previous logs, cvs and png files deleted!"


echo "===== Running final cell-type annotation ..................."

SUBsubDIR="dir_name"

## rm previous custom outputs
# rm ${CODEDIR}/${SUBDIR}/logs/06_jaccard_annotated_*.txt
# rm ${PROCESSEDIR}/${SUBDIR}/${SUBsubDIR}/Silhouette_*.cvs
# rm ${PLOTDIR}/${SUBDIR}/${SUBsubDIR}/*.pdf
# rm ${PLOTDIR}/${SUBDIR}/${SUBsubDIR}/*.png

echo "Previous outputs deleted!"


# ## Make Multiome UMAPS Broad and Fine
# 17_wnn_clustering_final_ct.sh
# 
# ## check ... 
# 18_wnn_profile_plots_final_ct.R
# 
# 19_geneExpr_plots_final_ct.R
# 
# 20_geneExpr_mean_ratio_final_ct.R
# 
# ## make heatmaps
# 21_heatmaps_markers.sh

# # Run dependency job
# id3=$(sbatch --parsable 05_rename_idents.sh)
# # Alternative command if you need to re-run the clustering
# # id3=$(sbatch --parsable --dependency=afterok:$id2 05_rename_idents.sh)
# echo $id3
# echo "Renaming Idents Done!!"

echo "**** Job ends ****"
date


} > $log_path 2>&1



## This script was made using slurmjobs version 1.2.5
## available from http://research.libd.org/slurmjobs/