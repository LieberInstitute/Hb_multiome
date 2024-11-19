#!/bin/bash
#SBATCH --partition=katun	
#SBATCH --job-name=run_all_outlier_detection
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=30GB						                                    # each job from the array will get its own private 30G to work with)
#SBATCH -o logs/run_all_outlier_detection.sh.%a.txt
#SBATCH -e logs/run_all_outlier_detection.sh.%a.txt
# SBATCH --mail-type=ALL

### Note CSC. SLURM Script in progress 

echo "**** Job starts ****"
echo "Run outlier detection in CellRangerARC-renalyze datasets"
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

########## (1-A) From SCRATCH (raw-counts) filter cells that PASS Outliers on `Cell Ranger ARC reanalyze` data, and save barcodes filtered

# ### Processed-data: ~processed-data/01_preprocessing_QC/cellrangerARC_reanalyze
# 
# ## remove RDS and CVS files with barcodes that passed outliers detected only on GEX assay
# # ~/01_preprocessing_QC/cellrangerARC_reanalyze
# "*_droplet_scores.rds"
# # ~/01_preprocessing_QC/cellrangerARC_reanalyze/csv_files
# "*_total_TrueCells.csv"
# "*_droplet_FDRcutoff.csv"
# "*_bc_PASS_isOutliers.csv"

# ########## Plots: ~/plots/01_preprocessing_QC/cellrangerARC_reanalyze
# 
# "*_droplet_qc.png" # Note. only generated when run from scratch
# "*_isOutliers_metrics.png" 
# 
# ## Run script
# #~/code/01_preprocessing_QC
# 
# 05_sce_seurat_emptydroplet.sh


########## (1-B) Filter cells that PASS Outliers directatly on `Cell Ranger ARC reanalyze` data, and save barcodes filtered

## Delete the logs/old-results, and re-submit EmptyDrops
cd ${CODEDIR}/01_preprocessing_QC

## logs
rm logs/*.txt

## proccessed dir
## CVS files with barcodes that passed outliers detected only on GEX assay
#rm ${PROCESSEDIR}/01_preprocessing_QC/cellrangerARC_reanalyze/csv_files/*_bc_PASS_isOutliers.csv

## Plots
rm ${PLOTDIR}/01_preprocessing_QC/cellrangerARC_reanalyze/plots_by_sample/*_isOutliers_per_sample.png
rm ${PLOTDIR}/01_preprocessing_QC/cellrangerARC_reanalyze/*.png

## Run script
#~/code/01_preprocessing_QC

sbatch 06_qc_scater_scran_metrics.sh
sbatch 06b_qc_scater_scran_metrics_integrated_plot.sh


########## (2) Filter cells that PASS Outliers directatly on `Cell Ranger ARC reanalyze` data, and save barcodes filtered

# # Hb_multiome/code/04_DiffExpr_Clustering_seurat/cellrangerARC_reanalyze_outliers
# # rm 
# "*_subset_Outliers.rds"
# "/cvs_files_markers/seurat.norm_counts_Harmony_All_subset_Outliers_cluster_info.csv"
# 
# 06_prepare_seurat_reanalyze_with_outliers.sh
# 
# 

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/