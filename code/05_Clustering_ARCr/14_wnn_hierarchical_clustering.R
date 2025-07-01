########################################################################
## Compute and Plot gene expression plots for top marker genes for one cell type 
## I used Deconvobuddies::findMarkers_1vAll(), a convenient wrapped (Scran/Deconvobuddies) to compute test.type="binom" (1vsALL)
## - Used: Default direction = "up".  Impacts p-values: if "up" genes with logFC < 0 will have p.value = 1
##
## Authors. CSC
## Date. Jun 30, 2025
##
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## DeconvoBuddies 1.1+ is not available for Bioconductor 3.18 (conda_R/4.3.x), so I am using conda_R/4.4.x instead
########################################################################

library("SingleCellExperiment")
#library("DeconvoBuddies")
library("purrr")
library("dplyr")
library("stringr")
library("ggplot2")
library("here")


## directories
inputSCE_Dir <- here(
    "processed-data", 
    "08_spatial_registration_vs_multiome_snRNA-seq"
)
processedDir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "14_wnn_hierarchical_clustering"
)
plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "14_wnn_hierarchical_clustering"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(processedDir)) {
    dir.create(processedDir)
}

#===============================================================================

## SingleCellExperiment object derived from Seurat rna modality
sce_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v4.rds"
sce_name <- here(inputSCE_Dir, sce_base_name)
sce_name
title_name <- str_extract(sce_base_name, regex("C\\.\\w*\\_r2"))
title_name

## Load SCE derived from Seurat WNN
sce <- readRDS(sce_name)

message("====== Verifications ======")
sce
# class: SingleCellExperiment 
# dim: 29690 55702 
# ...

## Verification
head(rownames(sce))
colnames(colData(sce))

## check cell types and counts
cluster_counts_df <- as.data.frame(table(sce$cluster_ann))
head(cluster_counts_df)

# check rowData ensembl names and gene id(s)
head(rowData(sce)$gene_symbol)
head(rowData(sce)$gene_id)

celltype_counts <- table(colData(sce)$cluster_ann)
as.data.frame(celltype_counts)


#===============================================================================



