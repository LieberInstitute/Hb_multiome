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
library("purrr")
library("dplyr")
library("stringr")
library("ggplot2")
library("here")

## directories

## clusters renamed for Spatial-Registration on Visium project

inputSCE_Dir <- "~/Habenula_Visium/processed-data/05_snRNA-seq_model_stats/"
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

source(here("code", "05_Clustering_ARCr", "get_mean_ratio_sparse.R"))
#get_mean_ratio_sparse

#===============================================================================

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

#===============================================================================
##  remove cell types with fewer than 10 cells
## Louise might be interested into add this filter into her get_mean_ratio() 
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
#Removing cell types with less<10 cells: C.42.no-match

##==============================================================================
## Get the mean ratio for each gene for each cell type defined in `cluster_ann`

# Debugging warning triggered by `get_mean_ratio()` 
# Warning message:
#     In asMethod(object) :
#     sparse->dense coercion: allocating vector of size 12.3 GiB

# verify "logcounts" is already a sparse matrix. OK
#class(assay(sce, "logcounts"))
# [1] "dgCMatrix"
# attr(,"package")
# [1] "Matrix"
#DeconvoBuddies:::get_mean_ratio  # See internal logic

# To avoid coercing a massive matrix, restrict the calculation to a subset of
# 10k genes most highly expressed genes
# keep_genes <- head(order(Matrix::rowMeans(assay(sce, "logcounts")), decreasing = TRUE), 2000) 
# sce_subset <- sce[keep_genes, ]
# sce_subset
# 
# ## support allocate a vector less than size 1.2 GiB
# marker_stats_ratio <- get_mean_ratio(
#     sce_subset,
#     assay_name = "logcounts",
#     cellType_col = "cluster_ann", 
#     gene_name = "gene_symbol",
#     gene_ensembl = "gene_id"
# )
# head(marker_stats_ratio, n=5)
# # verification
# top50_per_celltype <- marker_stats_ratio |>
#     filter(str_detect(cellType.target, "LHb|MHb")) |>
#     group_by(cellType.target) |>
#     slice_max(order_by = MeanRatio, n = 50, with_ties = FALSE) |>
#     ungroup()
# head(top50_per_celltype, n=10)
# 
# c("GPR151", "POUF4", "TAC3") %in% top50_lhb_mhb$gene
# top50_lhb_mhb %>%
#     filter(gene == "GPR151")

##==============================================================================
## inspect data
#head(marker_stats_ratio)
# # A tibble: 6 × 10
# gene    cellType.target     mean.target cellType.2nd        mean.2nd MeanRatio
# <chr>   <chr>                     <dbl> <chr>                  <dbl>     <dbl>
# 1 DLGAP2  C.25.undeterminated        1.36 C.04.undeterminated     1.13     1.19 
# 2 MALAT1  C.25.undeterminated        5.37 C.33.DD_LHb             5.67     0.946
# 3 MEG3    C.25.undeterminated        3.12 C.32.undeterminated     3.38     0.923
# 4 SNHG14  C.25.undeterminated        4.06 C.31.DD_Exit.Thal       4.42     0.919
# 5 SIPA1L1 C.25.undeterminated        1.06 C.31.DD_Exit.Thal       1.15     0.917
# 6 FTX     C.25.undeterminated        2.68 C.33.DD_LHb             3.14     0.856

##==============================================================================
## MeanRatio ranking is context-dependent on the gene universe you're testing
## - change the denominator in a relative ranking
## - MeanRatio = mean_in_target_celltype / max_mean_in_non_target_celltypes

##### Using 6000 genes 
# A tibble: 6 × 10
# gene    cellType.target mean.target cellType.2nd      mean.2nd MeanRatio
# <chr>   <chr>                 <dbl> <chr>                <dbl>     <dbl>
# 1 CALN1   C.05.DD_LHb           1.83  C.31.DD_Exit.Thal    1.33       1.38
# 2 TMTC4   C.05.DD_LHb           0.758 C.07.DD_MHb          0.577      1.31
# 3 CBLN2   C.05.DD_LHb           1.23  C.33.DD_LHb          0.972      1.27
# 4 SRGAP1  C.05.DD_LHb           2.05  C.10.DD_MHb          1.65       1.24
# 5 PRR16   C.05.DD_LHb           2.01  C.07.DD_MHb          1.66       1.21
# 6 L3MBTL4 C.05.DD_LHb           1.04  C.07.DD_MHb          0.866      1.20

##### Using 2000 genes 
# A tibble: 10 × 10
# gene    cellType.target mean.target cellType.2nd       mean.2nd MeanRatio
# <chr>   <chr>                 <dbl> <chr>                 <dbl>     <dbl>
# 1 CALN1   C.05.DD_LHb            1.83 C.31.DD_Exit.Thal      1.33      1.38
# 2 SRGAP1  C.05.DD_LHb            2.05 C.10.DD_MHb            1.65      1.24
# 3 PRR16   C.05.DD_LHb            2.01 C.07.DD_MHb            1.66      1.21
# 4 COL25A1 C.05.DD_LHb            3.31 C.16.DD_MHb            2.77      1.20
# 5 PBX3    C.05.DD_LHb            1.46 C.38.DD_Inhib.Thal     1.29      1.14
# 6 TMTC1   C.05.DD_LHb            2.33 C.07.DD_MHb            2.09      1.12


##==============================================================================
## Track MeanRatio Rank across gene subsets for Habenula genes. The goals are: 
# - Get an stable high MeanRatio and low rank across subset sizes. Indicates robust marker performance
# - If rank drops significantly when increasing feature space, that marker may be less specific

genes_of_interest <- c("GPR151", "TAC3", "POU4F1")

# Wrapper function to compute and extract ranks
get_marker_ranks <- function(sce, gene_subset, n, celltype_regex = "LHb|MHb",
                             genes_of_interest = NULL) {
    
    sce_sub <- sce[gene_subset, ]
    
    ratio_df <- get_mean_ratio_sparse(
        sce_sub,
        assay_name = "logcounts",
        cellType_col = "cluster_ann",
        gene_name = "gene_symbol",
        gene_ensembl = "gene_id"
    ) |>
        filter(str_detect(cellType.target, celltype_regex)) |>
        group_by(cellType.target) |>
        mutate(Rank = rank(-MeanRatio, ties.method = "first")) |>
        ungroup() |>
        mutate(SubsetSize = n)
    # Apply gene filter only if genes_of_interest is provided
    if (!is.null(genes_of_interest)) {
        ratio_df <- ratio_df |> filter(gene %in% genes_of_interest)
    }   
    
    return(ratio_df)
}

# Run across different gene subset sizes
sizes <- c(2000, 4000, 6000, 8000, 10000, 12000)

rank_results <- lapply(sizes, function(n) {
    gene_subset <- head(order(Matrix::rowMeans(assay(sce, "logcounts")), decreasing = TRUE), n)
    get_marker_ranks(sce, gene_subset, n, genes_of_interest = genes_of_interest)
})
length(rank_results)
rank_summary <- bind_rows(rank_results)

## Check how many rows were returned 
rank_summary$SubsetSize <- as.numeric(rank_summary$SubsetSize)
table(rank_summary$SubsetSize)
# 4000  6000  8000 10000 12000 
# 6    11    11    12    12 

## check clusters present
unique(rank_summary$cellType.target)
length(unique(rank_summary$cellType.target))
# [1] "C.11.DD_MHb" "C.05.DD_LHb" "C.07.DD_MHb" "C.10.DD_MHb" "C.30.DD_LHb"
# [6] "C.16.DD_MHb"
# only those 6 clusters had at least one of the genes_of_interest (GPR151, TAC3, POU4F1) ranked among the top markers

f_name <- here(plotDir ,"Track_MeanRatioRank_across_Hb_genes.pdf")
pdf(file = f_name, width = 12, height = 8)  # standard letter size

ggplot(rank_summary, aes(x = SubsetSize, y = MeanRatio, color = gene)) +
    geom_line(aes(group = interaction(gene, cellType.target))) +
    geom_point() +
    scale_x_continuous(breaks = sizes) +
    theme_minimal() +
    facet_wrap(~ cellType.target) +
    labs(title = "MeanRatio of Habenula marker genes across gene subset sizes",
         y = "MeanRatio", x = "Number of genes used")

dev.off()
# X-axis (SubsetSize): number of genes used in the input to get_mean_ratio_sparse()
# y-axis: ratio of average expression in the target cell type vs the most similar (highest mean) non-target cell type for that gene


##==============================================================================
## Get best-performing markers by stability (low variability in MeanRatio)

rank_results <- lapply(sizes, function(n) {
    gene_subset <- head(order(Matrix::rowMeans(assay(sce, "logcounts")), decreasing = TRUE), n)
    get_marker_ranks(sce, gene_subset, n)
})
length(rank_results)
names(rank_results) <- paste0("size_", sizes)
rank_summary <- bind_rows(rank_results)

## check clusters present
unique(rank_summary$cellType.target)
# [1] "C.11.DD_MHb" "C.23.DD_LHb" "C.05.DD_LHb" "C.18.DD_LHb" "C.14.DD_MHb"
# [6] "C.24.DD_LHb" "C.07.DD_MHb" "C.36.DD_MHb" "C.10.DD_MHb" "C.33.DD_LHb"
# [11] "C.30.DD_LHb" "C.16.DD_MHb" "C.40.DD_LHb"
length(unique(rank_summary$cellType.target))
# [1] 13

## Compute Marker Stability Metrics
stability_summary <- rank_summary |>
    group_by(gene, cellType.target) |>
    summarise(
        mean_ratio_mean = mean(MeanRatio, na.rm = TRUE),
        mean_ratio_min = min(MeanRatio, na.rm = TRUE),
        mean_ratio_max = max(MeanRatio, na.rm = TRUE),
        stability_range = mean_ratio_max - mean_ratio_min,
        mean_ratio_sd = sd(MeanRatio, na.rm = TRUE),
        .groups = "drop"
    )
head(stability_summary)
# gene  cellType.target mean_ratio_mean mean_ratio_min mean_ratio_max
# <chr> <chr>                     <dbl>          <dbl>          <dbl>
# 1 AAK1  C.05.DD_LHb               0.721          0.721          0.721
# 2 AAK1  C.07.DD_MHb               0.511          0.511          0.511
# 3 AAK1  C.10.DD_MHb               0.562          0.562          0.562
# 4 AAK1  C.11.DD_MHb               0.589          0.589          0.589
# 5 AAK1  C.14.DD_MHb               0.583          0.583          0.583
# 6 AAK1  C.16.DD_MHb               0.543          0.543          0.543

# Top 5 most stable (smallest range) per cluster
top_stable_markers <- stability_summary |>
    group_by(cellType.target) |>
    slice_min(order_by = stability_range, n = 3, with_ties = FALSE) |>
    ungroup()
head(top_stable_markers)

## Filter the Original Data for These Markers
rank_top_stable <- rank_summary |>
    semi_join(top_stable_markers, by = c("gene", "cellType.target"))

## Plot MeanRatio Over Subset Size

f_name <- here(plotDir ,"Top3_Stable_MarkerGenes_per_Cluster.pdf")
pdf(file = f_name, width = 12, height = 8)  # standard letter size

ggplot(rank_top_stable, aes(x = SubsetSize, y = MeanRatio, color = gene)) +
    geom_line(aes(group = interaction(gene, cellType.target))) +
    geom_point() +
    facet_wrap(~ cellType.target) +
    scale_x_continuous(breaks = unique(rank_summary$SubsetSize)) +
    theme_minimal() +
    labs(title = "Top 3 Stable Marker Genes per Cluster",
         y = "MeanRatio", x = "Number of Genes Used")

dev.off()


##==============================================================================

## save marker stats / rank_results
marker_ranks_6000 <- rank_results[["size_6000"]]
marker_ranks_8000 <- rank_results[["size_8000"]]
marker_ranks_10000 <- rank_results[["size_10000"]]
marker_ranks_12000 <- rank_results[["size_12000"]]
f_name <- here(processedDir, "marker_ranks_6k_12k.RData")
save(marker_ranks_6000, marker_ranks_8000, marker_ranks_10000, marker_ranks_12000, file = f_name)

message("Mean ratio results saved !!!")

##==============================================================================
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
# [1] "C.05.DD_LHb"         "C.07.DD_MHb"         "C.10.DD_MHb"        
# [4] "C.11.DD_MHb"         "C.14.DD_MHb"         "C.16.DD_MHb"        
# [7] "C.18.DD_LHb"         "C.23.DD_LHb"         "C.24.DD_LHb"        
# [10] "C.30.DD_LHb"         "C.33.DD_LHb"         "C.36.DD_MHb"        
# [13] "C.40.DD_LHb"         "C.01.undeterminated" "C.02.DD_Oligo"      
# [16] "C.03.undeterminated" "C.04.undeterminated" "C.06.DD_Exit.Thal"  
# [19] "C.08.undeterminated" "C.09.undeterminated" "C.12.undeterminated"
# [22] "C.13.no-match"       "C.15.DD_Exit.Thal"   "C.17.DD_Exit.Thal"  
# [25] "C.19.DD_Inhib.Thal"  "C.20.DD_Astrocyte"   "C.21.DD_Astrocyte"  
# [28] "C.22.undeterminated" "C.25.undeterminated" "C.26.DD_OPC"        
# [31] "C.27.DD_Microglia"   "C.28.DD_Inhib.Thal"  "C.29.DD_Endo"       
# [34] "C.31.DD_Exit.Thal"   "C.32.undeterminated" "C.34.DD_Oligo"      
# [37] "C.35.undeterminated" "C.37.undeterminated" "C.38.DD_Inhib.Thal" 
# [40] "C.39.DD_Inhib.Thal"  "C.41.DD_Microglia"  

# set as factor of sce[[cellType_col]] to match sorted_ct on the plots
sce$cluster_ann <- factor(sce$cluster_ann, levels = sorted_ct)

f_name <- here(plotDir ,"VPlot_mean_ratio_all_wnn_cluster_genes10k.pdf")
pdf(file = f_name, width = 8.5, height = 11)  # standard letter size

marker_stats <- rank_results[["size_10000"]]
head(marker_stats)
# A tibble: 6 × 12
# gene   cellType.target mean.target cellType.2nd        mean.2nd MeanRatio
# <chr>  <chr>                 <dbl> <chr>                  <dbl>     <dbl>
# 1 COX17  C.11.DD_MHb           0.688 C.30.DD_LHb            0.591      1.16
# 2 CHRNA3 C.11.DD_MHb           1.12  C.07.DD_MHb            0.965      1.16
# 3 NCS1   C.11.DD_MHb           0.858 C.03.undeterminated    0.802      1.07
# 4 GNB1   C.11.DD_MHb           1.30  C.14.DD_MHb            1.23       1.06
# 5 ITM2C  C.11.DD_MHb           1.74  C.10.DD_MHb            1.65       1.06
# 6 NRN1   C.11.DD_MHb           1.20  C.10.DD_MHb            1.13       1.06


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
        rank_col = "MeanRatio.rank",
        anno_col = "MeanRatio.anno",
        color_pal = NULL,
        plot_points = FALSE,
        ncol = 2
    )
    print(p1)
}

dev.off()

message("Mean ratio Plots for genes subset 10k saved !!!")



## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
