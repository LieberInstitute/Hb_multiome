#!/bin/bash
#SBATCH --partition=katun	
#SBATCH --job-name=Hb_reanalyze_all
#SBATCH -c 1
#SBATCH -t 1-00:00:00
#SBATCH --mem=50GB						                                    # each job from the array will get its own private 30G to work with)
#SBATCH --array=1-10                                           	  # change the array range acordingly with the rows listed in the array_targets.txt
#SBATCH --output=logs/Hb_reanalyze_all_%A_%a.out
#SBATCH --error=logs/Hb_reanalyze_all_%A_%a.out
# SBATCH --mail-type=ALL

echo "**** Job starts ****"
echo "Run on CellRanger-ARC reanalyze data sets"
date

echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job name (script): ${SLURM_JOB_NAME}"
echo "Hostname (computer node): ${HOSTNAME}"
echo ""
echo "Job id: ${SLURM_JOBID}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "Array/sample: $id"

## Relevant NOTE: Seurat objects should keep under the same conda version to avaid unexpected data conversion

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

######## Workflow for `Multiome 10x CellRanger ARC reanalyze scRNAseq + scATACseq data sets`

## Build Seurats and Plot basic QCs ########

cd ${CODEDIR}/01_preprocessing_QC

echo "Measure Quality Controls for GEX and ATAC assays from CellRanger-ARC data using Seurat & Signac packages"
echo "/01_preprocessing_QC/01_preprocessing_GEX_ATAC.sh"
sbatch ${CODEDIR}/01_preprocessing_QC/01_preprocessing_GEX_ATAC.sh


echo "Build Seurat and calculate initial Quality Controls from `CellRanger-ARC reanalyze`"
echo "/01_preprocessing_QC02_preprocessing_GEX_ATAC_reanalyze.sh"
## remove old files
rm logs/*.out
#rm ${PROCESSEDIR}/01_preprocessing_QC/cellrangerARC_reanalyze/*_reanalysis.rds
rm ${PROCESSEDIR}/01_preprocessing_QC/cellrangerARC_reanalyze/*metrics_summary.csv
rm ${PLOTDIR}/01_preprocessing_QC/cellrangerARC_reanalyze/Seurat_Signac_QC_metrics/*.png
sbatch 02_preprocessing_GEX_ATAC_reanalyze.sh # ?
sbatch 04_preprocessing_QCs_reanalyze.R # check difference btw these scripts


## Combine Seurats   ########################

# ## Delete the logs/old-results, and re-submit seurat builder
# cd ${CODEDIR}/02_merge_seurats
# # mkdir -p logs ## Create the logs directory if it doesn't exist
# rm logs/01_merge_seurats_job_loop_*.txt
# rm ${PROCESSEDIR}/02_merge_seurats/*seurat.combined.data_counts.rds
# rm ${PROCESSEDIR}/02_merge_seurats/*seurat.combined.norm_counts.rds
# ## These are dependency jobs-arrays 
# #sbatch 01_merge_seurats_job_loop.sh
# id1=$(sbatch --parsable 01_merge_seurats_job_loop.sh)
# 
# 
# ######## Perform CCA and Harmony corrections  ########
# 
# ## Delete the logs/old-results, and re-submit seurat builder
# #cd ${CODEDIR}/02_merge_seurats
# # mkdir -p logs ## Create the logs directory if it doesn't exist
# rm logs/02_integrate_seurats_job_loop*.txt
# rm ${PROCESSEDIR}/02_merge_seurats/*CCA.rds
# rm ${PROCESSEDIR}/02_merge_seurats/*Harmony.rds
# rm ${PROCESSEDIR}/02_merge_seurats/*PCA.rds
# ## These are dependency jobs-arrays from id1 
# id2=$(sbatch --parsable --dependency=afterok:$id1 02_integrate_seurats_job_loop.sh)
# 
# 
# ######## Integrate some plots of `RNA` side from multiome  ########
# 
# ## Delete the logs/old-results, and re-submit seurat builder
# #cd ${CODEDIR}/02_merge_seurats
# # mkdir -p logs ## Create the logs directory if it doesn't exist
# rm logs/03_merge_seurat_gex_plts*.txt
# rm ${PROCESSEDIR}/02_merge_seurats/*seurat.combined.norm_counts_all_plots_GEX.pdf
# ## These are dependency job from id1
# id2=$(sbatch --parsable --dependency=afterok:$id1 03_merge_seurat_gex_plts.sh)
# 
# 
# ## Future steps here
# ## I am working on the sh script of this chunk yet
# 
# 
# ######## Calculate DEG in batch corrected data and pseudobulk datasets  ########
# 
# rm logs/aggregateExpression_genes*.txt
# rm ${PROCESSEDIR}/03_pseudobulking/*_pseudobulk.rds
# rm ${PROCESSEDIR}/03_pseudobulking/*_subset.rds
# rm ${PROCESSEDIR}/03_pseudobulking/cvs_files_markers/*_Allmarkers_min*cells.csv
# ## These is a dependency job from id2. It requires `Harmony or CCA` Seurat RDS datasets
# id3=$(sbatch --parsable --dependency=afterok:$id2 01_aggregateExpression_genes.sh)
# #01_aggregateExpression_genes.R
# 
# 
# ######## Read Seurat clusters to search cell-types based on a custom marker gene list   ########
# 04_DiffExpr_Clustering_seurat/01_Hb_celltypes_from_seurat_obj.R
# 
# ######## Calculate and summarize percentage statistics of pre-selected cell types in different categories   ########
# ## You need to update clusters of interes before run this script
# # allT <- c(unique(df_mdT[["seurat_clusters"]]))
# # hb <- c(1, 7, 8, 13, 18, 19) ## Habenula clusters
# # neu <- c(1, 5, 6, 7, 8, 11, 13, 18, 19) ## Neuron clusters
# 04_DiffExpr_Clustering_seurat/02_cell_types_percentages.R
# 
# 
# echo "**** Job ends ****"
# date

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/

