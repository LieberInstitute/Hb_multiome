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
library("dplyr")
library("stringr")
library("here")

## directories

## clusters renamed for Spatial-Registration on Visium project

inputSCE_Dir <- "~/Habenula_Visium/processed-data/05_snRNA-seq_model_stats/"
outputCSV_Dir <- here(
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
if (!dir.exists(outputCSV_Dir)) {
    dir.create(outputCSV_Dir)
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
## briefly check cell types and counts
cluster_counts_df <- as.data.frame(table(sce$cluster_ann))
head(cluster_counts_df)
#                   Var1 Freq
# 1 C.01.undeterminated 3906
# 2       C.02.DD_Oligo 3371
# 3 C.03.undeterminated 3271
# 4 C.04.undeterminated 2911
# 5         C.05.DD_LHb 2771
# 6   C.06.DD_Exit.Thal 2667

# check rowData ensembl names and gene id(s)
head(rowData(sce)$gene_symbol)
head(rowData(sce)$gene_id)

##  remove cell types with fewer than 10 cells
celltypes <- colData(sce)$cluster_ann
celltype_counts <- table(celltypes)
low_ct <- names(celltype_counts[celltype_counts <= 10])

## Louise might be interested into add this filter into the package
if (length(low_ct)==TRUE) {
    message("Removing cell types with less<10 cells: ", low_ct)
    # Remove cell types with <= 10 cells
    valid_types <- names(celltype_counts[celltype_counts > 10])
    sce <- sce[, colData(sce)$cluster_ann %in% valid_types]
}

##==============================================================================
## Get the mean ratio for each gene for each cell type defined in `cluster_ann`

# Debugging warning triggered by `get_mean_ratio()` 
# Warning message:
#     In asMethod(object) :
#     sparse->dense coercion: allocating vector of size 12.3 GiB

# verify "logcounts" is already a sparse matrix. OK
class(assay(sce, "logcounts"))
# [1] "dgCMatrix"
# attr(,"package")
# [1] "Matrix"

# To avoid coercing a massive matrix, restrict the calculation to a subset of 
# 3000 genes most highly expressed genes
keep_genes <- head(order(Matrix::rowMeans(assay(sce, "logcounts")), decreasing = TRUE), 10000)
sce_subset <- sce[keep_genes, ]

marker_stats <- get_mean_ratio(
    sce_subset,
    assay_name = "logcounts",
    cellType_col = "cluster_ann", 
    gene_name = "gene_symbol",
    gene_ensembl = "gene_id"
)

##==============================================================================


## inspect data
head(marker_stats)
# # A tibble: 6 × 10
# gene    cellType.target     mean.target cellType.2nd        mean.2nd MeanRatio
# <chr>   <chr>                     <dbl> <chr>                  <dbl>     <dbl>
# 1 DLGAP2  C.25.undeterminated    1.36 C.04.undeterminated     1.13     1.19 
# 2 MALAT1  C.25.undeterminated        5.37 C.33.DD_LHb             5.67     0.946
# 3 MEG3    C.25.undeterminated        3.12 C.32.undeterminated     3.38     0.923
# 4 SNHG14  C.25.undeterminated        4.06 C.31.DD_Exit.Thal       4.42     0.919
# 5 SIPA1L1 C.25.undeterminated        1.06 C.31.DD_Exit.Thal       1.15     0.917
# 6 FTX     C.25.undeterminated        2.68 C.33.DD_LHb             3.14     0.856

filtered_marker_stats <- marker_stats |>
    filter(cellType.target == "C.05.DD_LHb")
filtered_marker_stats
# gene    cellType.target mean.target cellType.2nd      mean.2nd MeanRatio
# <chr>   <chr>                 <dbl> <chr>                <dbl>     <dbl>
#     1 MSC-AS1 C.05.DD_LHb           0.971 C.24.DD_LHb          0.396      2.45
# 2 GALR1   C.05.DD_LHb           1.00  C.36.DD_MHb          0.667      1.50
# 3 CALN1   C.05.DD_LHb           1.83  C.31.DD_Exit.Thal    1.33       1.38
# 4 TMTC4   C.05.DD_LHb           0.758 C.07.DD_MHb          0.577      1.31
# 5 CBLN2   C.05.DD_LHb           1.23  C.33.DD_LHb          0.972      1.27
# 6 SRGAP1  C.05.DD_LHb           2.05  C.10.DD_MHb          1.65       1.24
# 7 PRR16   C.05.DD_LHb           2.01  C.07.DD_MHb          1.66       1.21
# 8 L3MBTL4 C.05.DD_LHb           1.04  C.07.DD_MHb          0.866      1.20
# 9 COL25A1 C.05.DD_LHb           3.31  C.16.DD_MHb          2.77       1.20
# 10 ZNF235  C.05.DD_LHb           0.529 C.33.DD_LHb          0.458      1.16

## save marker stats
marker_stats
save(marker_stats, file = here(outputCSV_Dir, sprintf("marker_stats_MeanRatio_%s.Rdata", cluster)))

## plots the top n marker genes for a specified cell type based off of the stats table from get_mean_ratio()

# Prepare data
# sort unique_ct so that clusters containing "MH" or "LH" appear first
# valid_types = after remove ct<10 cells
sorted_ct <- valid_types|>
    sort() |>
    tibble(cluster = _) |>
    mutate(priority = str_detect(cluster, "MHb|LHb")) |>
    arrange(desc(priority), cluster) |>
    pull(cluster)
sorted_ct

f_name <- here(plotDir ,"VPlot_mean_ratio_all_wnn_cluster.pdf")
pdf(file = f_name, width = 8.5, height = 11)  # standard letter size

# Loop through each habenula cluster
for (ct in sorted_ct) {
    message("Plotting wnn cluster: ", ct)
    p1 <- plot_marker_express(
        sce,
        stats = marker_stats,
        cellType_col = "cluster_ann",
        cell_type = ct,
        gene_col = "gene",
        n_genes = 10,
        # rank_col = "MeanRatio.rank",
        # anno_col = "MeanRatio.anno",
        color_pal = NULL,
        plot_points = FALSE,
        ncol = 2
    )
    print(p1)
}

dev.off()

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
