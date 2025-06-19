########################################################################
## Plot gene expression violin plots for top marker genes for one cell type 
## Use Mean-Ratio
##
## Authors. CSC
## Date. Jun 19, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
##
# All Seurat objects were created using module load conda_R/4.3.x to preserve chromatin integrity.
# If you're not working with the ATAC modality and encounter installation issues with other packages
## -—such as DeconvoBuddies, which is not available for Bioconductor version '3.18', you may consider using conda_R/4.4.x instead
# Relevant Note: 
# - Avoid updating your Seurat objects under this version, as it may compromise chromatin integrity.
# - Warning message: package ‘DeconvoBuddies’ is not available for Bioconductor version '3.18'
########################################################################

library("Seurat")
library("DeconvoBuddies")
library("here")

## directories

## clusters renamed for Spatial-Registration on Visium project

inputRDS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "05_rename_idents"
)
inputCVS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "02_Hb_celltypes_from_seurat_reanalyze_v3",
    "cvs_files_markers"
)
plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "08_wnn_geneExp_plts_annotated"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}

## Load input with RDS wnn to compare
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+"))
# C.leiden_lsi_r2_renamed_visium

seurat_name <- here(inputRDS_Dir, Seurat_base_name)
title_name <- str_extract(seurat_name, regex("C\\.\\w*\\_r2"))

# Load Seurat
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
# verification
levels(SeuratOBJ)
colnames(SeuratOBJ@meta.data)

## make slim Seurat with ony RNA modality
SeuratOBJ <- DietSeurat(SeuratOBJ, 
                        assays = "RNA",
                        Reductions = NULL)
SeuratOBJ
head(rownames(SeuratOBJ[["RNA"]]))

## convert Seurat object into sce
sce <- as.SingleCellExperiment(SeuratOBJ)
## verification
sce

## find marker genes with get_mean_ratio
marker_stats <- get_mean_ratio(
    sce_DLPFC_example,
    cellType_col = "cellType_broad_hc",
    gene_name = "gene_name",
    gene_ensembl = "gene_id"
)



## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
