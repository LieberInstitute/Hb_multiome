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
# Relevant Notes: 
# (1) Avoid updating your Seurat objects under this version, as it may compromise chromatin integrity.
#     - Warning message: package ‘DeconvoBuddies’ is not available for Bioconductor version '3.18'
#     - ERROR: this R is version 4.3.2, package 'DeconvoBuddies' requires R >=  4.4.0
# (2) Alternatively, install development version from GitHub 
#.    - BiocManager::install("LieberInstitute/DeconvoBuddies")
########################################################################

library("SingleCellExperiment")
library("DeconvoBuddies") # [1] ‘0.99.39’
library("stringr")
library("here")

## directories

## clusters renamed for Spatial-Registration on Visium project

inputSCE_Dir <- "~/Habenula_Visium/processed-data/05_snRNA-seq_model_stats/"
# inputCVS_Dir <- here(
#     "processed-data",
#     "05_Clustering_ARCr",
#     "02_Hb_celltypes_from_seurat_reanalyze_v3",
#     "cvs_files_markers"
# )
plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "13_wnn_geneExp_plt_mean_ratio_annotated"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}

## set hard path to `sce` object derived from Seurat multimodal dataset
sce_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v3.rds"
sce_name <- paste0(inputSCE_Dir, sce_base_name)
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

## nuclei are classified in to cell types
table(sce$cluster_ann)
# C.01.undeterminated       C.02.DD_Oligo C.03.undeterminated C.04.undeterminated 
# 3906                3371                3271                2911 
# C.05.DD_LHb   C.06.DD_Exit.Thal         C.07.DD_MHb C.08.undeterminated 
# 2771                2667                2610                2537 
# C.09.undeterminated         C.10.DD_MHb         C.11.DD_MHb C.12.undeterminated 
# 2430                2344                2212                2187 
# C.13.no-match         C.14.DD_MHb   C.15.DD_Exit.Thal         C.16.DD_MHb 
# 2036                2026                1625                1607 
# C.17.DD_Exit.Thal         C.18.DD_LHb  C.19.DD_Inhib.Thal   C.20.DD_Astrocyte 
# 1587                1535                1380                1343 
# C.21.DD_Astrocyte C.22.undeterminated         C.23.DD_LHb         C.24.DD_LHb 
# 1341                1327                1269                 825 
# C.25.undeterminated         C.26.DD_OPC   C.27.DD_Microglia  C.28.DD_Inhib.Thal 
# 707                 638                 587                 543 
# C.29.DD_Endo         C.30.DD_LHb   C.31.DD_Exit.Thal C.32.undeterminated 
# 343                 213                 209                 196 
# C.33.DD_LHb       C.34.DD_Oligo C.35.undeterminated         C.36.DD_MHb 
# 186                 184                 165                 145 
# C.37.undeterminated  C.38.DD_Inhib.Thal  C.39.DD_Inhib.Thal         C.40.DD_LHb 
# 111                 105                  90                  84 
# C.41.DD_Microglia       C.42.no-match 
# 76                   2 

# check rowData ensembl names and gene id(s)
head(rowData(sce)$gene_symbol)
head(rowData(sce)$gene_id)

##  remove cell types with fewer than 10 cells
celltypes <- colData(sce)$cluster_ann
celltype_counts <- table(celltypes)
low_ct <- names(celltype_counts[celltype_counts <= 10])

if (length(low_ct)==TRUE) {
    message("Removing cell types with less<10 cells: ", low_ct)
    # Remove cell types with <= 10 cells
    valid_types <- names(celltype_counts[celltype_counts > 10])
    sce <- sce[, colData(sce)$cluster_ann %in% valid_types]
}

## Get the mean ratio for each gene for each cell type defined in `cluster_ann`
# specify rowData col names for gene_name and gene_ensembl
marker_stats <- get_mean_ratio(
    sce,
    cellType_col = "cluster_ann", 
    gene_name = "gene_symbol",
    gene_ensembl = "gene_id"
)

head(marker_stats)
# # A tibble: 6 × 8
# gene    cellType.target     mean.target cellType.2nd        mean.2nd MeanRatio
# <chr>   <chr>                     <dbl> <chr>                  <dbl>     <dbl>
# 1 DLGAP2  C.25.undeterminated        1.36 C.04.undeterminated     1.13     1.19 
# 2 MALAT1  C.25.undeterminated        5.37 C.33.DD_LHb             5.67     0.946
# 3 MEG3    C.25.undeterminated        3.12 C.32.undeterminated     3.38     0.923
# 4 SNHG14  C.25.undeterminated        4.06 C.31.DD_Exit.Thal       4.42     0.919
# 5 SIPA1L1 C.25.undeterminated        1.06 C.31.DD_Exit.Thal       1.15     0.917
# 6 FTX     C.25.undeterminated        2.68 C.33.DD_LHb             3.14     0.856


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
