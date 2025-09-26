########################################################################
## Call LinkPeaks() on pseudobulk at Mid cell-type level
## Parse Seurat by ct due memory issues
## 
## Authors. CSC
## Date. Sep 07, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("BSgenome.Hsapiens.UCSC.hg38")
library("GenomicRanges") 
library("purrr")
library("tidyverse")
library("tidyr")
library("stringr")
library("here")


## ATAC function's helper used globally
source(here("code", "06_peak_calling", "multiome_custom_functions", "multiome_idents_normalization_helper.R"))
# ls()

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

p_met = "spearman"
w_size_num <- 5e5 
w_size = "5e5" # for filenames only 
resolution_level = "Mid" 
#cluster_name <- "Endo"

## read input arguments
args = commandArgs(trailingOnly = TRUE)
cluster_name <- args[2]
# cluster_name = "Inhib.Thal"
if (is.na(cluster_name) || !nzchar(cluster_name)) stop("Missing cluster_name argument")

## Check/create directories
inputRDS_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "12_pseudobulk_MACS2"
)
output_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "13_pseudobulk_LinkPeaks_MACS2_split_ct",
    "links_ct_merged"
)

if (!dir.exists(output_Dir)) {
    dir.create(output_Dir)
}


##==============================================================================
## Load Seurat / macs peaks / filtered genes. And make verification

SeuratOBJ <- readRDS(here(inputRDS_Dir, "Mid_pseudobulk.spearman.5e5.rds"))
pb_obj <- SeuratOBJ
pb_obj
# An object of class Seurat 
# 371023 features across 169 samples within 2 assays 
# Active assay: ATAC_macs2_pseudo (355127 features, 337546 variable features)
# 2 layers present: counts, data
# 1 other assay present: RNA
# 1 dimensional reduction calculated: lsi

DefaultAssay(pb_obj) <- "RNA"
levels(pb_obj)
# [1] "Astrocyte"  "Endo"       "Excit.Thal" "Inhib.Thal" "LHb.1"     
# [6] "LHb.1.3"    "LHb.1.3.4"  "LHb.2.7"    "LHb.4"      "LHb.7"     
# [11] "MHb.1"      "MHb.1.2"    "MHb.2"      "MHb.3"      "Microglia" 
# [16] "Oligo"      "OPC"        "Thal" 

message("Assays present: ", paste(Assays(pb_obj), collapse=", "))
# Assays present: RNA, ATAC_macs2_pseudo


##==============================================================================
# genome for RegionStats
genome <- BSgenome.Hsapiens.UCSC.hg38

# derive cluster labels
if (!"cluster_pb" %in% colnames(pb_obj@meta.data)) {
    pb_obj$cluster_pb <- sub("[-_].*$", "", colnames(pb_obj))
}
pb_obj$cluster_pb <- trimws(pb_obj$cluster_pb)
Idents(pb_obj) <- pb_obj$cluster_pb
clusters_pb <- levels(Idents(pb_obj))

message(length(clusters_pb), " Clusters in pseudobulk: ", paste(clusters_pb, collapse=", "))

stopifnot(cluster_name %in% clusters_pb)

# test all genes in pb_obj
genes_for_lp <- rownames(pb_obj[["RNA"]])


##=====================================================================

# Use the pseudobulk ATAC assay that exists
atac_assay_name <- if ("ATAC_macs2_pseudo" %in% Assays(pb_obj)) "ATAC_macs2_pseudo" else stop("ATAC_macs2_pseudo does not exist!")
DefaultAssay(pb_obj) <- atac_assay_name

# peaks from assay and peak->cluster map (from metadata)
peaks_gr <- granges(pb_obj[[atac_assay_name]])
head(peaks_gr)
# GRanges object with 6 ranges and 1 metadata column:
#     seqnames        ranges strand |         peak_called_in
# <Rle>     <IRanges>  <Rle> |            <character>
# [1]     chr1 181329-181534      * |             Inhib.Thal
# [2]     chr1 191217-191619      * | OPC,Oligo,Inhib.Thal..
# [3]     chr1 629146-629354      * |              Astrocyte
# [4]     chr1 629811-630032      * | LHb.7,Astrocyte,Olig..
# [5]     chr1 630189-630389      * |             Oligo,Endo
# [6]     chr1 632189-632410      * |              Astrocyte
peaks_gr$peak_id <- Signac::GRangesToString(peaks_gr)

# MACS2 GRanges stored as metadata column 'peak_called_in' in peaks_gr
stopifnot("peak_called_in" %in% colnames(mcols(peaks_gr)))

peaks_map <- as.data.frame(peaks_gr) |>
    transmute(peak_id = peak_id,
              peak_called_in = as.character(peak_called_in)) |>
    mutate(peak_called_in = str_split(peak_called_in, "\\s*,\\s*")) |>
    tidyr::unnest(peak_called_in) |>
    mutate(peak_called_in = trimws(peak_called_in)) |>
    filter(!is.na(peak_called_in), peak_called_in != "") |>
    distinct()

head(peaks_map)
# # A tibble: 6 × 2
# peak_id            peak_called_in
# <chr>              <chr>         
# 1 chr1-181329-181534 Inhib.Thal    
# 2 chr1-191217-191619 OPC           
# 3 chr1-191217-191619 Oligo         
# 4 chr1-191217-191619 Inhib.Thal    
# 5 chr1-191217-191619 LHb.4         
# 6 chr1-191217-191619 MHb.2 

peaks_by_cluster <- function(cluster) {
    peaks_map |> filter(peak_called_in == cluster) |> pull(peak_id) |> unique()
}

# subset pseudobulk to the requested cluster
seurat_subset <- subset(pb_obj, idents = cluster_name)
n_samples <- ncol(seurat_subset)
message("Processing ", cluster_name, " (", n_samples, " pseudobulk samples)")

if (n_samples < 3) stop("Too few samples for cluster ", cluster_name, "")


##==============================================================================

# keep only this cluster’s MACS2-called peaks present in assay
cluster_peaks <- intersect(peaks_by_cluster(cluster_name),
                           rownames(seurat_subset[[atac_assay_name]]))
if (length(cluster_peaks) == 0) stop("No MACS2 peaks for ", cluster_name, " in assay.")

# drop ultra-sparse peaks in this cluster (improves stability)
counts_mat <- GetAssayData(
    seurat_subset, 
    assay = atac_assay_name, 
    layer = "counts")[cluster_peaks, , drop = FALSE]

min_cells_sub <- max(3, floor(0.05 * n_samples))   # ≥5% or at least 3
keep_peaks_sub <- rownames(counts_mat)[Matrix::rowSums(counts_mat > 0) >= min_cells_sub]
if (length(keep_peaks_sub) == 0) stop("No peaks pass support filter in ", cluster_name, ".")


# Rebuild assay with kept peaks (preserves annotation cleanly)
peak_ranges <- Signac::StringToGRanges(keep_peaks_sub)

new_assay  <- CreateChromatinAssay(
    counts = counts_mat[keep_peaks_sub, , drop = FALSE],
    ranges = peak_ranges,
    annotation = tryCatch(Annotation(pb_obj[[atac_assay_name]]), error = function(e) NULL)
)
# ChromatinAssay data with 5225 features for 10 cells
# Variable features: 0 
# Genome: 
#     Annotation present: TRUE 
# Motifs present: FALSE 
# Fragment files: 0 

seurat_subset[[atac_assay_name]] <- new_assay
# Typically when assigning a new ChromatinAssay, a warning tells the new assay doesn’t have exactly the same feature set as the old one (ATAC_macs2_pseudo). That’s expected after merged or filtered peaks. Harmless.

DefaultAssay(seurat_subset) <- atac_assay_name

# RegionStats (attach GC%, width) — needed for bias correction
seurat_subset <- RegionStats(
    object = seurat_subset,
    assay  = atac_assay_name,
    genome = genome
)
# occasionally the peaks set contains seqnames not present in the genome package (e.g., random contigs, chrM if filtered, unplaced scaffolds). Those peaks simply won’t get stats assigned. Harmless

seurat_subset <- global_rebuild_atac_normalization(seurat_subset, atac_assay_name)
# Seurat auto-detects that only 9 SVD are valid (effective rank of the pb count matrix)
# silently trims to 9 as I have only ~10 pseudobulk groups, and I can’t extract 50 components anyway.

# Genes present in this subset’s RNA
genes_sub <- intersect(genes_for_lp, rownames(seurat_subset[["RNA"]]))
if (length(genes_sub) == 0) stop("No RNA genes to test in subset.")

message("Summary:",
        "\n  Cell-type = ", cluster_name,
        "\n  Peaks = ", length(keep_peaks_sub),
        "\n  Genes = ", length(genes_sub),
        "\n  min.cells = ", min_cells_sub,
        "\n  distance = ", w_size_num)
seurat_subset
# Summary:
#     Cell-type = LHb.1.3
# Peaks = 3409
# Genes = 15896
# min.cells = 3
# distance = 5e+05
# An object of class Seurat 

## ============================================================================/

message("Computing peak-gene correlations on pseudobulk by cell-type ")

DefaultAssay(seurat_subset) <- "RNA"  # not strictly required but conventional

# retain all links for downstream multiple testing correction (FDR/HB)
seurat_subset <- LinkPeaks(
    object           = seurat_subset,
    peak.assay       = atac_assay_name,
    expression.assay = "RNA",
    genes.use        = genes_sub,
    distance         = w_size_num,
    min.cells        = min_cells_sub,
    pvalue_cutoff    = 1,   # keep everything / keep only links with pvalue <= pvalue.cutoff
    score_cutoff     = 0,   # keep both positive and negative scores
    method           = p_met
)

# Extract links and compute FDR
lk_gr <- Links(seurat_subset[[atac_assay_name]])
lk_df <- as.data.frame(lk_gr)

# Add FDR column
lk_df$FDR <- p.adjust(lk_df$pvalue, method = "BH")
lk_df$cluster <- cluster_name

message("Links found: ", nrow(lk_df))
# Links found: 28251
if (nrow(lk_df) > 0) {
    print(head(lk_df[, c("peak", "gene", "score", "pvalue", "FDR")], 5))
}

# Save as CSV
out_csv <- here(output_Dir, paste0(resolution_level, "_", cluster_name,
                                   "_pseudobulk_link_peak_genes.csv"))
write.csv(lk_df, out_csv, row.names = FALSE)
message("Links: ", nrow(lk_df), " (saved: ", basename(out_csv), ")")

# Save subset if desired
f_name <- paste0(resolution_level, "_", cluster_name,
                 "_pseudobulk_seurat_subset.rds")
output_Dir <- here(output_Dir, "seurats_rds", f_name)
saveRDS(seurat_subset, file = output_Dir)

message("All done!!!")


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


# > ## Reproducibility information
# > library("sessioninfo")
# > print("Reproducibility information:")
# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2025-09-26 10:59:38 EDT"
# > proc.time()
# user  system elapsed 
# 13.743   1.414 649.591 
# > options(width = 120)
# > session_info()
# 0 (R 4.4.2)
# bitops                        1.0-9     2024-10-03 [2] CRAN (R 4.4.1)
# BSgenome                    * 1.74.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# BSgenome.Hsapiens.UCSC.hg38 * 1.4.5     2024-05-29 [2] Bioconductor
# cli                           3.6.5     2025-04-23 [2] CRAN (R 4.4.3)
# cluster                       2.1.8     2024-12-11 [3] CRAN (R 4.4.3)
# codetools                     0.2-20    2024-03-31 [3] CRAN (R 4.4.3)
# colorspace                    2.1-1     2024-07-26 [2] CRAN (R 4.4.1)
# cowplot                       1.1.3     2024-01-22 [2] CRAN (R 4.4.0)
# crayon                        1.5.3     2024-06-20 [2] CRAN (R 4.4.1)
# curl                          6.3.0     2025-06-06 [1] CRAN (R 4.4.3)
# data.table                    1.17.6    2025-06-17 [1] CRAN (R 4.4.3)
# DelayedArray                  0.32.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# deldir                        2.0-4     2024-02-28 [2] CRAN (R 4.4.0)
# dichromat                     2.0-0.1   2022-05-02 [2] CRAN (R 4.4.0)
# digest                        0.6.37    2024-08-19 [2] CRAN (R 4.4.1)
# dotCall64                     1.2       2024-10-04 [2] CRAN (R 4.4.1)
# dplyr                       * 1.1.4     2023-11-17 [2] CRAN (R 4.4.0)
# farver                        2.1.2     2024-05-13 [2] CRAN (R 4.4.0)
# fastDummies                   1.7.5     2025-01-20 [2] CRAN (R 4.4.2)
# fastmap                       1.2.0     2024-05-15 [2] CRAN (R 4.4.0)
# fastmatch                     1.1-6     2024-12-23 [2] CRAN (R 4.4.2)
# fitdistrplus                  1.2-2     2025-01-07 [2] CRAN (R 4.4.2)
# forcats                     * 1.0.0     2023-01-29 [2] CRAN (R 4.4.0)
# future                        1.49.0    2025-05-09 [2] CRAN (R 4.4.3)
# future.apply                  1.11.3    2024-10-27 [2] CRAN (R 4.4.2)
# generics                      0.1.4     2025-05-09 [2] CRAN (R 4.4.3)
# GenomeInfoDb                * 1.42.3    2025-01-27 [2] Bioconductor 3.20 (R 4.4.2)
# GenomeInfoDbData              1.2.13    2024-10-01 [2] Bioconductor
# GenomicAlignments             1.42.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# GenomicRanges               * 1.58.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# ggplot2                     * 3.5.2     2025-04-09 [2] CRAN (R 4.4.3)
# ggrepel                       0.9.6     2024-09-07 [2] CRAN (R 4.4.1)
# ggridges                      0.5.6     2024-01-23 [2] CRAN (R 4.4.0)
# ggvenn                      * 0.1.10    2023-03-31 [1] CRAN (R 4.4.3)
# globals                       0.18.0    2025-05-08 [2] CRAN (R 4.4.3)
# glue                          1.8.0     2024-09-30 [2] CRAN (R 4.4.1)
# goftest                       1.2-3     2021-10-07 [2] CRAN (R 4.4.0)
# gridExtra                     2.3       2017-09-09 [2] CRAN (R 4.4.0)
# gtable                        0.3.6     2024-10-25 [2] CRAN (R 4.4.2)
# here                        * 1.0.1     2020-12-13 [2] CRAN (R 4.4.0)
# hms                           1.1.3     2023-03-21 [2] CRAN (R 4.4.0)
# htmltools                     0.5.8.1   2024-04-04 [2] CRAN (R 4.4.0)
# htmlwidgets                   1.6.4     2023-12-06 [2] CRAN (R 4.4.0)
# httpuv                        1.6.16    2025-04-16 [2] CRAN (R 4.4.3)
# httr                          1.4.7     2023-08-15 [2] CRAN (R 4.4.0)
# ica                           1.0-3     2022-07-08 [2] CRAN (R 4.4.0)
# igraph                        2.1.4     2025-01-23 [2] CRAN (R 4.4.2)
# IRanges                     * 2.40.1    2024-12-05 [2] Bioconductor 3.20 (R 4.4.2)
# irlba                         2.3.5.1   2022-10-03 [2] CRAN (R 4.4.0)
# jsonlite                      2.0.0     2025-03-27 [2] CRAN (R 4.4.3)
# KernSmooth                    2.23-26   2025-01-01 [3] CRAN (R 4.4.3)
# later                         1.4.2     2025-04-08 [2] CRAN (R 4.4.3)
# lattice                       0.22-6    2024-03-20 [3] CRAN (R 4.4.3)
# lazyeval                      0.2.2     2019-03-15 [2] CRAN (R 4.4.0)
# lifecycle                     1.0.4     2023-11-07 [2] CRAN (R 4.4.0)
# listenv                       0.9.1     2024-01-29 [2] CRAN (R 4.4.0)
# lmtest                        0.9-40    2022-03-21 [2] CRAN (R 4.4.0)
# lubridate                   * 1.9.4     2024-12-08 [2] CRAN (R 4.4.2)
# magrittr                      2.0.3     2022-03-30 [2] CRAN (R 4.4.0)
# MASS                          7.3-65    2025-02-28 [3] CRAN (R 4.4.3)
# Matrix                        1.7-2     2025-01-23 [3] CRAN (R 4.4.3)
# MatrixGenerics                1.18.1    2025-01-09 [2] Bioconductor 3.20 (R 4.4.2)
# matrixStats                   1.5.0     2025-01-07 [2] CRAN (R 4.4.2)
# mime                          0.13      2025-03-17 [2] CRAN (R 4.4.3)
# miniUI                        0.1.2     2025-04-17 [2] CRAN (R 4.4.3)
# nlme                          3.1-167   2025-01-27 [3] CRAN (R 4.4.3)
# parallelly                    1.44.0    2025-05-07 [2] CRAN (R 4.4.3)
# patchwork                     1.3.0     2024-09-16 [2] CRAN (R 4.4.1)
# pbapply                       1.7-2     2023-06-27 [2] CRAN (R 4.4.0)
# pillar                        1.10.2    2025-04-05 [2] CRAN (R 4.4.3)
# pkgconfig                     2.0.3     2019-09-22 [2] CRAN (R 4.4.0)
# plotly                        4.10.4    2024-01-13 [2] CRAN (R 4.4.0)
# plyr                          1.8.9     2023-10-02 [2] CRAN (R 4.4.0)
# png                           0.1-8     2022-11-29 [2] CRAN (R 4.4.0)
# polyclip                      1.10-7    2024-07-23 [2] CRAN (R 4.4.1)
# progressr                     0.15.1    2024-11-22 [2] CRAN (R 4.4.2)
# promises                      1.3.3     2025-05-29 [1] CRAN (R 4.4.3)
# purrr                       * 1.0.4     2025-02-05 [2] CRAN (R 4.4.2)
# R6                            2.6.1     2025-02-15 [2] CRAN (R 4.4.2)
# RANN                          2.6.2     2024-08-25 [2] CRAN (R 4.4.1)
# RColorBrewer                  1.1-3     2022-04-03 [2] CRAN (R 4.4.0)
# Rcpp                          1.0.14    2025-01-12 [2] CRAN (R 4.4.2)
# RcppAnnoy                     0.0.22    2024-01-23 [2] CRAN (R 4.4.0)
# RcppHNSW                      0.6.0     2024-02-04 [2] CRAN (R 4.4.0)
# RcppRoll                      0.3.1     2024-07-07 [1] CRAN (R 4.4.2)
# RCurl                         1.98-1.17 2025-03-22 [2] CRAN (R 4.4.3)
# readr                       * 2.1.5     2024-01-10 [2] CRAN (R 4.4.0)
# reshape2                      1.4.4     2020-04-09 [2] CRAN (R 4.4.0)
# restfulr                      0.0.15    2022-06-16 [2] CRAN (R 4.4.0)
# reticulate                    1.42.0    2025-03-25 [2] CRAN (R 4.4.3)
# rjson                         0.2.23    2024-09-16 [2] CRAN (R 4.4.1)
# rlang                         1.1.6     2025-04-11 [2] CRAN (R 4.4.3)
# ROCR                          1.0-11    2020-05-02 [2] CRAN (R 4.4.0)
# rprojroot                     2.0.4     2023-11-05 [2] CRAN (R 4.4.0)
# Rsamtools                     2.22.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# RSpectra                      0.16-2    2024-07-18 [2] CRAN (R 4.4.1)
# rtracklayer                 * 1.66.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# Rtsne                         0.17      2023-12-07 [2] CRAN (R 4.4.0)
# S4Arrays                      1.6.0     2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# S4Vectors                   * 0.44.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# scales                        1.4.0     2025-04-24 [2] CRAN (R 4.4.3)
# scattermore                   1.2       2023-06-12 [2] CRAN (R 4.4.0)
# sctransform                   0.4.2     2025-04-30 [2] CRAN (R 4.4.3)
# sessioninfo                 * 1.2.3     2025-02-05 [2] CRAN (R 4.4.2)
# Seurat                      * 5.3.0     2025-04-23 [2] CRAN (R 4.4.3)
# SeuratObject                * 5.1.0     2025-04-22 [2] CRAN (R 4.4.3)
# shiny                         1.10.0    2024-12-14 [2] CRAN (R 4.4.2)
# Signac                      * 1.14.0    2024-08-21 [1] CRAN (R 4.4.2)
# sp                          * 2.2-0     2025-02-01 [2] CRAN (R 4.4.2)
# spam                          2.11-1    2025-01-20 [2] CRAN (R 4.4.2)
# SparseArray                   1.6.2     2025-02-20 [2] Bioconductor 3.20 (R 4.4.3)
# spatstat.data                 3.1-6     2025-03-17 [2] CRAN (R 4.4.3)
# spatstat.explore              3.4-2     2025-03-21 [2] CRAN (R 4.4.3)
# spatstat.geom                 3.4-1     2025-05-20 [2] CRAN (R 4.4.3)
# spatstat.random               3.4-1     2025-05-20 [2] CRAN (R 4.4.3)
# spatstat.sparse               3.1-0     2024-06-21 [2] CRAN (R 4.4.1)
# spatstat.univar               3.1-3     2025-05-08 [2] CRAN (R 4.4.3)
# spatstat.utils                3.1-4     2025-05-15 [2] CRAN (R 4.4.3)
# stringi                       1.8.7     2025-03-27 [2] CRAN (R 4.4.3)
# stringr                     * 1.5.1     2023-11-14 [2] CRAN (R 4.4.0)
# SummarizedExperiment          1.36.0    2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# survival                      3.8-3     2024-12-17 [3] CRAN (R 4.4.3)
# tensor                        1.5       2012-05-05 [2] CRAN (R 4.4.0)
# tibble                      * 3.3.0     2025-06-08 [1] CRAN (R 4.4.3)
# tidyr                       * 1.3.1     2024-01-24 [2] CRAN (R 4.4.0)
# tidyselect                    1.2.1     2024-03-11 [2] CRAN (R 4.4.0)
# tidyverse                   * 2.0.0     2023-02-22 [2] CRAN (R 4.4.0)
# timechange                    0.3.0     2024-01-18 [2] CRAN (R 4.4.0)
# tzdb                          0.5.0     2025-03-15 [2] CRAN (R 4.4.3)
# UCSC.utils                    1.2.0     2024-10-29 [2] Bioconductor 3.20 (R 4.4.2)
# uwot                          0.2.3     2025-02-24 [2] CRAN (R 4.4.3)
# vctrs                         0.6.5     2023-12-01 [2] CRAN (R 4.4.0)
# viridisLite                   0.4.2     2023-05-02 [2] CRAN (R 4.4.0)

