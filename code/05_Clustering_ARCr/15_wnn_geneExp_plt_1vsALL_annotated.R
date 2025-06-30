########################################################################
## Compute and Plot gene expression plots for top marker genes for one cell type 
## I used Deconvobuddies::findMarkers_1vAll(), a convenient wrapped (Scran/Deconvobuddies) to compute test.type="binom" (1vsALL)
## 
## Authors. CSC
## Date. Jun 30, 2025
##
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## DeconvoBuddies 1.1+ is not available for Bioconductor 3.18 (conda_R/4.3.x), so I am using conda_R/4.4.x instead
########################################################################

library("SingleCellExperiment")
library("DeconvoBuddies")
library("purrr")
library("dplyr")
library("stringr")
library("ggplot2")
library("here")


## directories

#inputSCE_Dir <- "~/Habenula_Visium/processed-data/05_snRNA-seq_model_stats/"
inputSCE_Dir <- here(
    "processed-data", 
    "08_spatial_registration_vs_multiome_snRNA-seq"
)
processedDir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "13_wnn_geneExp_plt_mean_ratio_annotated"
)
plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "13_wnn_geneExp_plt_mean_ratio_annotated"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(processedDir)) {
    dir.create(processedDir)
}

# no longer required - fixed on recent deconvobuddies release
#source(here("code", "05_Clustering_ARCr", "get_mean_ratio_sparse.R"))
#get_mean_ratio_sparse

#===============================================================================

## SingleCellExperiment object derived from Seurat rna modality
sce_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v4.rds"
sce_name <- here(inputSCE_Dir, sce_base_name)
sce_name
title_name <- str_extract(sce_base_name, regex("C\\.\\w*\\_r2"))
title_name

# Load Seurat
sce <- readRDS(sce_name)
## verification
sce
# class: SingleCellExperiment 
# dim: 29690 55702 
# metadata(0):
#     assays(3): counts logcounts scaledata
# rownames(29690): MIR1302-2HG FAM138A ... AC007325.4 AC007325.2
# rowData names(2): gene_symbol gene_id
# colnames(55702): S04_AAACAGCCAGAATGAC-1 S04_AAACAGCCAGCAAGGC-1 ...
# S09_TTTGTTGGTCATGCAA-1 S09_TTTGTTGGTTGTTCAC-1
# colData names(31): orig.ident nCount_RNA ... cluster_ann ident
# reducedDimNames(0):
#     mainExpName: RNA
# altExpNames(0):

head(rownames(sce))

colnames(colData(sce))

## briefly check cell types and counts
cluster_counts_df <- as.data.frame(table(sce$cluster_ann))
head(cluster_counts_df)

# check rowData ensembl names and gene id(s)
head(rowData(sce)$gene_symbol)
head(rowData(sce)$gene_id)

#===============================================================================
## Prepare data: remove cell types with fewer than 10 cells

celltypes <- colData(sce)$cluster_ann
celltype_counts <- table(celltypes)
low_ct <- names(celltype_counts[celltype_counts <= 10])

# filter into the get_mean_ratio() function
if (length(low_ct)==TRUE) {
    message("Removing cell types with less<10 cells: ", low_ct)
    # Remove cell types with <= 10 cells
    valid_types <- names(celltype_counts[celltype_counts > 10])
    sce <- sce[, colData(sce)$cluster_ann %in% valid_types]
}


