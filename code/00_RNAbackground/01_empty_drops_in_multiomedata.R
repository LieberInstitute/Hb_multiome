## Code used as reference
## https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/03_droplet_qc.R
# Adapted by HT. LIEBER `Habenula_Pilot` project (code/04_snRNA-seq/01_get_droplet_scores.R)
# https://github.com/LieberInstitute/Habenula_Pilot/blob/fa452f307a32ab063417d8ab8002917f82c1f703/code/04_snRNA-seq/01_get_droplet_scores.R

## run command in slurm:
## $ sbatch slurm_emptydrops_V2.sh
## INPUT with samples names to process is on ../01_create_array_target.txt


########################################################################
## Runs empty drops in multimome cellranger-ARC data.
## Goal: Distinguish between droplets containing cells and ambient RNA in a droplet-based single-cell RNA sequencing experiment
## .     Get preliminary statistics
##      Processes non-empty cells
## .     Create and save a table with the statistics
## .     Save the single-cell object processed
## .     Create and save a knee plot
## Authors. CSC
## Date. August 14th, 2023
## Last.md Jan 19th, 2024
########################################################################


library(Seurat)
library(SingleCellExperiment)
library(DropletUtils) # Functions for scRNA-seq data from droplet technologies such as 10X Genomics
library(dplyr)
library(ggplot2)
library(here)
library(sessioninfo)

here::here()
set.seed(100)


############  Initials ############

# Call functions to create and handle meta-data to seurat objects
source(here("code/functions_custom", "remote_file_caller.R"))

FDR_cutoff <- 0.001

## commandArgs scans the arguments which have been supplied when the current R script was invoked (from shell sh)
sample_tmp <- commandArgs(trailingOnly = TRUE)
# sample_tmp <- args[1]
# for testing: sample_tmp <- '5S_Hb_KDM,human'

sample_data <- unlist(strsplit(sample_tmp, ","))
s_sample <- sample_data[[1]]

message("Processing sample: ", s_sample)

## Break to avoid processing old samples
if (s_sample %in% c("S1_Hb_KDM", "S1_Hb_KDM")) {
    print("Skipe old sample. ")
    stop()
}

## Read the raw_feature_bc_matrix.h5
sample_path <- get_raw_barcode_mtx(s_sample)
message("Reading data from ", sample_path) # ../cellrangerARC/S1_Hb_KDM/outs/raw_feature_bc_matrix.h5
h5_raw_path <- Read10X_h5(sample_path) # dgCMatrix data. Barcodes for columns and genes by rows
# head(h5_raw_path, n=1)                                 # Sparse mtx has the 2 slots (gene expression and peaks)
# Extract the 'Gene Expression' matrix only
raw.sce <- h5_raw_path$`Gene Expression` # SingleCellExperiment data.
head(raw.sce, n = 1)

## Get total number of cells in the gene expression assay
totalCells <- length(Cells(raw.sce))

## Compute barcode rank statistics and identify the knee and inflection points on the total count curve
bcRanks <- barcodeRanks(raw.sce, fit.bounds = c(10, 1e3))
## Barcode rank range.
colnames(bcRanks)
range((bcRanks$rank))
range((bcRanks$total))
## Get knee value and add hundred points to perform a more stringent threshold.
knee_lower <- metadata(bcRanks)$knee + 100
## Get the inflection point. Signs change.
inflection <- metadata(bcRanks)$inflection

message("Sample: ", s_sample, "\nKnee lower value: ", knee_lower, "\nInflection point:", inflection)

## Run emptyDrops w/ knee + 100 to make it more strident
## emptyDrops distinguish between droplets containing cells and ambient RNA in a droplet-based single-cell RNA experiment
st <- Sys.time()
message(Sys.time(), " Starting emptyDrops")
sce.out <- DropletUtils::emptyDrops(
    raw.sce,
    niters = 30000, # number of iterations. It use the Monte Carlo p-value
    lower = knee_lower # numeric scalar specifying the lower bound on the total UMI count
)
en <- Sys.time() - st

message(paste0(" Processing time: ", en))
# head.matrix(sce.out,n=5)

# Get significant TRUE cells based on the FDR cutoff
cells_FALSE <- 0
cells_FT <- 0
cells_TRUE <- 0
signif_TRUE <- 0
##  add arbitrary margins on a multidimensional array
cells_stat <- addmargins(table(
    Signif = sce.out$FDR <= FDR_cutoff,
    Limited = sce.out$Limited,
    useNA = "ifany"
))
cells_stat

# get specific values from a confusion mtx
cells_FALSE <- cells_stat[1, 1]
cells_FT <- cells_stat[2, 1]
cells_TRUE <- cells_stat[2, 2]
signif_TRUE <- cells_stat[2, 4]

# Calculate non Emptydroplets value and percentage related
nonEmptydroplets <- (sce.out |> as.data.frame() |> filter(FDR < FDR_cutoff) |> summarise(n = n()))$n
per.nonemptydroplets <- ((nonEmptydroplets * 100) / totalCells) # , digits = 4)

# Build a table with the stats applied and the outputs gotten
# Table with the cellranger-arc `gene expression statistics`
tab_stats <- matrix(c(
    s_sample, totalCells, knee_lower, inflection, FDR_cutoff,
    cells_FALSE, cells_FT, cells_TRUE, signif_TRUE,
    nonEmptydroplets, per.nonemptydroplets
), ncol = 11, byrow = TRUE)
colnames(tab_stats) <- c(
    "Sample_name", "Total_cells", "Knee_lower", "Inflection", "FDR_cutoff",
    "Signif.FALSE", "F/T", "T/T", "Signif.TRUE",
    "Non.Emptydroplets", "%Emptydroplets"
)
tab_stats <- as.table(tab_stats)

# > tab_stats
# Sample_name Total_cells Knee_lower Inflection FDR_cutoff Signif.FALSE F/T
# A S1_Hb_KDM   720339      956        268        0.001      407          224
# T/T  Signif.TRUE Non.Emptydroplets %Emptydroplets
# A 7351 7575        7575              1.05158821055087

# Export the table to CSV for further analysis
message(paste0(" Process completed. Saving stats and plots"))

s_file_name <- here("processed-data", "00_RNAbackground", paste0(s_sample, "_empty_droplets_stats.csv"))
write.csv(tab_stats, file = s_file_name, quote = FALSE, row.names = FALSE)

# Droplet Elbow plot
s_file_name <- here("plots", "00_RNAbackground", paste0(s_sample, "_empty_droplets_knee_plot.png"))

define_theme <- function(size = 15) {
    theme_bw() +
        theme(text = element_text(size = size))
}
# Prepare data frame with additional FDR column
droplet_elbow_data <- as.data.frame(bcRanks) %>%
    mutate(FDR = sce.out$FDR)

# Define parameters
knee_meta <- metadata(bcRanks)$knee
knee_lower_label <- paste0("Knee est 'lower' (", knee_lower, ")")
second_knee_label <- paste0("Second Knee (", knee_meta, ")")
title <- paste0("Sample: ", s_sample)
subtitle <- nonEmptydroplets # n_cell_anno

# Create ggplot object
droplet_elbow_plot <- droplet_elbow_data %>%
    ggplot(aes(x = rank, y = total, color = FDR < FDR_cutoff)) +
    # Define points
    geom_point(alpha = 0.5, size = 1) +
    # Define lines and annotations
    geom_hline(yintercept = knee_meta, linetype = "dotted", color = "gray") +
    annotate("text", x = 10, y = knee_meta, label = second_knee_label, vjust = -1, color = "gray") +
    geom_hline(yintercept = knee_lower, linetype = "dashed") +
    annotate("text", x = 10, y = knee_lower, label = knee_lower_label, vjust = -0.5) +
    # Define scales
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    # Define labels
    labs(
        x = "Barcode Rank",
        y = "Total UMI Counts",
        title = title,
        subtitle = subtitle,
        color = paste("FDR <", FDR_cutoff)
    ) +
    define_theme() +
    theme(legend.position = "bottom")

ggsave(droplet_elbow_plot, filename = s_file_name)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


# [1] "Reproducibility information:"
# > Sys.time()
# [1] "2024-08-08 14:50:50 EDT"
# > proc.time()
# user   system  elapsed
# 1055.326   15.723 9628.191
# > options(width = 120)
# > session_info()
# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 Patched (2024-02-08 r85876)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-08-08
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3.x/bin/pandoc
#
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version    date (UTC) lib source
# abind                  1.4-5      2016-07-21 [2] CRAN (R 4.3.2)
# beachmat               2.18.0     2023-10-24 [2] Bioconductor
# Biobase              * 2.62.0     2023-10-24 [2] Bioconductor
# BiocGenerics         * 0.48.1     2023-11-01 [2] Bioconductor
# BiocParallel           1.36.0     2023-10-24 [2] Bioconductor
# bit                    4.0.5      2022-11-15 [2] CRAN (R 4.3.2)
# bit64                  4.0.5      2020-08-30 [2] CRAN (R 4.3.2)
# bitops                 1.0-7      2021-04-24 [2] CRAN (R 4.3.2)
# cli                    3.6.2      2023-12-11 [2] CRAN (R 4.3.2)
# cluster                2.1.6      2023-12-01 [3] CRAN (R 4.3.2)
# codetools              0.2-19     2023-02-01 [3] CRAN (R 4.3.2)
# colorspace             2.1-0      2023-01-23 [2] CRAN (R 4.3.2)
# cowplot                1.1.3      2024-01-22 [2] CRAN (R 4.3.2)
# crayon                 1.5.2      2022-09-29 [2] CRAN (R 4.3.2)
# data.table             1.15.0     2024-01-30 [2] CRAN (R 4.3.2)
# DelayedArray           0.28.0     2023-10-24 [2] Bioconductor
# DelayedMatrixStats     1.24.0     2023-10-24 [2] Bioconductor
# deldir                 2.0-2      2023-11-23 [2] CRAN (R 4.3.2)
# digest                 0.6.34     2024-01-11 [2] CRAN (R 4.3.2)
# dotCall64              1.1-1      2023-11-28 [2] CRAN (R 4.3.2)
# dplyr                * 1.1.4      2023-11-17 [2] CRAN (R 4.3.2)
# dqrng                  0.3.2      2023-11-29 [2] CRAN (R 4.3.2)
# DropletUtils         * 1.22.0     2023-10-24 [2] Bioconductor
# edgeR                  4.0.14     2024-01-29 [2] Bioconductor 3.18 (R 4.3.2)
# ellipsis               0.3.2      2021-04-29 [2] CRAN (R 4.3.2)
# fansi                  1.0.6      2023-12-08 [2] CRAN (R 4.3.2)
# farver                 2.1.1      2022-07-06 [2] CRAN (R 4.3.2)
# fastDummies            1.7.3      2023-07-06 [2] CRAN (R 4.3.2)
# fastmap                1.1.1      2023-02-24 [2] CRAN (R 4.3.2)
# fitdistrplus           1.1-11     2023-04-25 [2] CRAN (R 4.3.2)
# future                 1.33.1     2023-12-22 [2] CRAN (R 4.3.2)
# future.apply           1.11.1     2023-12-21 [2] CRAN (R 4.3.2)
# generics               0.1.3      2022-07-05 [2] CRAN (R 4.3.2)
# GenomeInfoDb         * 1.38.5     2023-12-28 [2] Bioconductor 3.18 (R 4.3.2)
# GenomeInfoDbData       1.2.11     2024-02-09 [2] Bioconductor
# GenomicRanges        * 1.54.1     2023-10-29 [2] Bioconductor
# ggplot2              * 3.4.4      2023-10-12 [2] CRAN (R 4.3.2)
# ggrepel                0.9.5      2024-01-10 [2] CRAN (R 4.3.2)
# ggridges               0.5.6      2024-01-23 [2] CRAN (R 4.3.2)
# globals                0.16.2     2022-11-21 [2] CRAN (R 4.3.2)
# glue                   1.7.0      2024-01-09 [2] CRAN (R 4.3.2)
# goftest                1.2-3      2021-10-07 [2] CRAN (R 4.3.2)
# gridExtra              2.3        2017-09-09 [2] CRAN (R 4.3.2)
# gtable                 0.3.4      2023-08-21 [2] CRAN (R 4.3.2)
# HDF5Array              1.30.0     2023-10-24 [2] Bioconductor
# hdf5r                  1.3.9      2024-01-14 [2] CRAN (R 4.3.2)
# here                 * 1.0.1      2020-12-13 [2] CRAN (R 4.3.2)
# htmltools              0.5.7      2023-11-03 [2] CRAN (R 4.3.2)
# htmlwidgets            1.6.4      2023-12-06 [2] CRAN (R 4.3.2)
# httpuv                 1.6.14     2024-01-26 [2] CRAN (R 4.3.2)
# httr                   1.4.7      2023-08-15 [2] CRAN (R 4.3.2)
# ica                    1.0-3      2022-07-08 [2] CRAN (R 4.3.2)
# igraph                 2.0.1.9008 2024-02-09 [2] Github (igraph/rigraph@39158c6)
# IRanges              * 2.36.0     2023-10-24 [2] Bioconductor
# irlba                  2.3.5.1    2022-10-03 [2] CRAN (R 4.3.2)
# jsonlite               1.8.8      2023-12-04 [2] CRAN (R 4.3.2)
# KernSmooth             2.23-22    2023-07-10 [3] CRAN (R 4.3.2)
# later                  1.3.2      2023-12-06 [2] CRAN (R 4.3.2)
# lattice                0.22-5     2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval               0.2.2      2019-03-15 [2] CRAN (R 4.3.2)
# leiden                 0.4.3.1    2023-11-17 [2] CRAN (R 4.3.2)
# lifecycle              1.0.4      2023-11-07 [2] CRAN (R 4.3.2)
# limma                  3.58.1     2023-10-31 [2] Bioconductor
# listenv                0.9.1      2024-01-29 [2] CRAN (R 4.3.2)
# lmtest                 0.9-40     2022-03-21 [2] CRAN (R 4.3.2)
# locfit                 1.5-9.8    2023-06-11 [2] CRAN (R 4.3.2)
# magrittr               2.0.3      2022-03-30 [2] CRAN (R 4.3.2)
# MASS                   7.3-60.0.1 2024-01-13 [3] CRAN (R 4.3.2)
# Matrix                 1.6-5      2024-01-11 [3] CRAN (R 4.3.2)
# MatrixGenerics       * 1.14.0     2023-10-24 [2] Bioconductor
# matrixStats          * 1.2.0      2023-12-11 [2] CRAN (R 4.3.2)
# mime                   0.12       2021-09-28 [2] CRAN (R 4.3.2)
# miniUI                 0.1.1.1    2018-05-18 [2] CRAN (R 4.3.2)
# munsell                0.5.0      2018-06-12 [2] CRAN (R 4.3.2)
# nlme                   3.1-164    2023-11-27 [3] CRAN (R 4.3.2)
# parallelly             1.36.0     2023-05-26 [2] CRAN (R 4.3.2)
# patchwork              1.2.0      2024-01-08 [2] CRAN (R 4.3.2)
# pbapply                1.7-2      2023-06-27 [2] CRAN (R 4.3.2)
# pillar                 1.9.0      2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig              2.0.3      2019-09-22 [2] CRAN (R 4.3.2)
# plotly                 4.10.4     2024-01-13 [2] CRAN (R 4.3.2)
# plyr                   1.8.9      2023-10-02 [2] CRAN (R 4.3.2)
# png                    0.1-8      2022-11-29 [2] CRAN (R 4.3.2)
# polyclip               1.10-6     2023-09-27 [2] CRAN (R 4.3.2)
# progressr              0.14.0     2023-08-10 [2] CRAN (R 4.3.2)
# promises               1.2.1      2023-08-10 [2] CRAN (R 4.3.2)
# purrr                  1.0.2      2023-08-10 [2] CRAN (R 4.3.2)
# R.methodsS3            1.8.2      2022-06-13 [2] CRAN (R 4.3.2)
# R.oo                   1.26.0     2024-01-24 [2] CRAN (R 4.3.2)
# R.utils                2.12.3     2023-11-18 [2] CRAN (R 4.3.2)
# R6                     2.5.1      2021-08-19 [2] CRAN (R 4.3.2)
# ragg                   1.2.7      2023-12-11 [2] CRAN (R 4.3.2)
# RANN                   2.6.1      2019-01-08 [2] CRAN (R 4.3.2)
# RColorBrewer           1.1-3      2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp                   1.0.12     2024-01-09 [2] CRAN (R 4.3.2)
# RcppAnnoy              0.0.22     2024-01-23 [2] CRAN (R 4.3.2)
# RcppHNSW               0.6.0      2024-02-04 [2] CRAN (R 4.3.2)
# RCurl                  1.98-1.14  2024-01-09 [2] CRAN (R 4.3.2)
# reshape2               1.4.4      2020-04-09 [2] CRAN (R 4.3.2)
# reticulate             1.35.0     2024-01-31 [2] CRAN (R 4.3.2)
# rhdf5                  2.46.1     2023-11-29 [2] Bioconductor 3.18 (R 4.3.2)
# rhdf5filters           1.14.1     2023-11-06 [2] Bioconductor
# Rhdf5lib               1.24.2     2024-02-07 [2] Bioconductor 3.18 (R 4.3.2)
# rlang                  1.1.3      2024-01-10 [2] CRAN (R 4.3.2)
# ROCR                   1.0-11     2020-05-02 [2] CRAN (R 4.3.2)
# rprojroot              2.0.4      2023-11-05 [2] CRAN (R 4.3.2)
# RSpectra               0.16-1     2022-04-24 [2] CRAN (R 4.3.2)
# Rtsne                  0.17       2023-12-07 [2] CRAN (R 4.3.2)
# S4Arrays               1.2.0      2023-10-24 [2] Bioconductor
# S4Vectors            * 0.40.2     2023-11-23 [2] Bioconductor 3.18 (R 4.3.2)
# scales                 1.3.0      2023-11-28 [2] CRAN (R 4.3.2)
# scattermore            1.2        2023-06-12 [2] CRAN (R 4.3.2)
# sctransform            0.4.1      2023-10-19 [2] CRAN (R 4.3.2)
# scuttle                1.12.0     2023-10-24 [2] Bioconductor
# sessioninfo          * 1.2.2      2021-12-06 [2] CRAN (R 4.3.2)
# Seurat               * 5.0.1      2023-11-17 [2] CRAN (R 4.3.2)
# SeuratObject         * 5.0.1      2023-11-17 [2] CRAN (R 4.3.2)
# shiny                  1.8.0      2023-11-17 [2] CRAN (R 4.3.2)
# SingleCellExperiment * 1.24.0     2023-10-24 [2] Bioconductor
# sp                   * 2.1-3      2024-01-30 [2] CRAN (R 4.3.2)
# spam                   2.10-0     2023-10-23 [2] CRAN (R 4.3.2)
# SparseArray            1.2.3      2023-12-25 [2] Bioconductor 3.18 (R 4.3.2)
# sparseMatrixStats      1.14.0     2023-10-24 [2] Bioconductor
# spatstat.data          3.0-4      2024-01-15 [2] CRAN (R 4.3.2)
# spatstat.explore       3.2-6      2024-02-01 [2] CRAN (R 4.3.2)
# spatstat.geom          3.2-8      2024-01-26 [2] CRAN (R 4.3.2)
# spatstat.random        3.2-2      2023-11-29 [2] CRAN (R 4.3.2)
# spatstat.sparse        3.0-3      2023-10-24 [2] CRAN (R 4.3.2)
# spatstat.utils         3.0-5      2024-06-17 [1] CRAN (R 4.3.2)
# statmod                1.5.0      2023-01-06 [2] CRAN (R 4.3.2)
# stringi                1.8.3      2023-12-11 [2] CRAN (R 4.3.2)
# stringr                1.5.1      2023-11-14 [2] CRAN (R 4.3.2)
# SummarizedExperiment * 1.32.0     2023-10-24 [2] Bioconductor
# survival               3.5-7      2023-08-14 [3] CRAN (R 4.3.2)
# systemfonts            1.0.5      2023-10-09 [2] CRAN (R 4.3.2)
# tensor                 1.5        2012-05-05 [2] CRAN (R 4.3.2)
# textshaping            0.3.7      2023-10-09 [2] CRAN (R 4.3.2)
# tibble                 3.2.1      2023-03-20 [2] CRAN (R 4.3.2)
# tidyr                  1.3.1      2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect             1.2.0      2022-10-10 [2] CRAN (R 4.3.2)
# utf8                   1.2.4      2023-10-22 [2] CRAN (R 4.3.2)
