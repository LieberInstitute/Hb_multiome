
### Note CSC. SLURM Script in progress 


########## (1-A) From SCRATCH (raw-counts) filter cells that PASS Outliers on `Cell Ranger ARC reanalyze` data, and save barcodes filtered

### Processed-data: ~processed-data/01_preprocessing_QC/cellrangerARC_reanalyze

## remove RDS and CVS files with barcodes that passed outliers detected only on GEX assay
# ~/01_preprocessing_QC/cellrangerARC_reanalyze
"*_droplet_scores.rds"
# ~/01_preprocessing_QC/cellrangerARC_reanalyze/csv_files
"*_total_TrueCells.csv"
"*_droplet_FDRcutoff.csv"
"*_bc_PASS_isOutliers.csv"

########## Plots: ~/plots/01_preprocessing_QC/cellrangerARC_reanalyze

"*_droplet_qc.png" # Note. only generated when run from scratch
"*_isOutliers_metrics.png" 

## Run script
#~/code/01_preprocessing_QC

05_sce_seurat_emptydroplet.sh
06_qc_scater_scran_metrics.sh 



########## (1-B) Filter cells that PASS Outliers directatly on `Cell Ranger ARC reanalyze` data, and save barcodes filtered

## remove CVS files with barcodes that passed outliers detected only on GEX assay
# ~/01_preprocessing_QC/cellrangerARC_reanalyze/csv_files
"*_bc_PASS_isOutliers.csv"

########## Plots: ~/plots/01_preprocessing_QC/cellrangerARC_reanalyze

"*_isOutliers_metrics.png" 

## Run script
#~/code/01_preprocessing_QC

06_qc_scater_scran_metrics.sh 



########## (2) Filter cells that PASS Outliers directatly on `Cell Ranger ARC reanalyze` data, and save barcodes filtered

# Hb_multiome/code/04_DiffExpr_Clustering_seurat
06_prepare_seurat_reanalyze_with_outliers.sh


