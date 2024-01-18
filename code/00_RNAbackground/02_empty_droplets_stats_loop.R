########################################################################
## Runs empty drops in cellranger-ARC data. 
## This script is designed for running as a job array in SGE or SLURM envs
## Goal:     
##.     Get preliminary statistics
##      Processes non-empty cells
##.     Create and save a table with the statistics
##.     Save the single-cell object processed
##.     Create and save a knee plot
##
## Implementation done in all the cellranger-ARC available datasets
##         
## Authors. CSC
## Date. August 14th, 2023
## Last.md August 17th, 2023
########################################################################

# Adapted from the LIEBER `Habenula_Pilot` project (code/04_snRNA-seq/01_get_droplet_scores.R) 
# https://github.com/LieberInstitute/Habenula_Pilot/blob/fa452f307a32ab063417d8ab8002917f82c1f703/code/04_snRNA-seq/01_get_droplet_scores.R
# qrsh -l mem_free=50G,h_vmem=50G

#### Load libraries ######
library(Seurat)
library(SingleCellExperiment)
library(DropletUtils)   # Functions for scRNA-seq data from droplet technologies such as 10X Genomics
library(dplyr)
library(ggplot2)
library(here)
library(sessioninfo)

here::here()

# Call functions to create and handle meta-data to seurat objects
source(here("code/functions_custom", "remote_file_caller.R"))

############  Initials ############
# The array target txt file is created with the 01_create_array_target.
# Run first if necessary calling: 
#     source(here("code", "01_create_array_tarjet.R"))    
# When written this script, the list of samples included in the array job are:
#     lst_samples <- list('hippo42_1', 'hippo42_4','1_HPC_KDM', '2_HPC_KDM', '3_HPC_KDM')

## commandArgs scans the arguments which have been supplied when the current R script was invoked (from shell sh)
sample_tmp <- commandArgs(trailingOnly = TRUE)
#sample_tmp <- args[1]
#sample_tmp <- '2_HPC_KDM,human'  # testing HUMAN tissue
#sample_tmp <- '3_HPC_KDM,mouse'  # testing MOUSE tissue 
sample_data = unlist(strsplit(sample_tmp,","))

s_sample <- sample_data[[1]]
#s_tissue <- sample_data[[2]]

# Read cellranger-ARC raw data from different paths (disks). 
sample_path <- get_raw_barcode_mtx(s_sample)
message("Reading data from ", sample_path)

#i <- 1              #control headers to create the table just in the first row
FDR_cutoff <- 0.001
set.seed(100)

# Accessing raw_data with h5 from Seurat objects takes longer time ...
st <- Sys.time()  # 17 seconds
h5_raw_path <- Read10X_h5(sample_path)                  # Get a dgCMatrix with bc for columns and genes by rows
#head(h5_raw_path, n=1)     #read just the gene expression slot
# Extract the 'Gene Expression' matrix from the h5 object
raw.sce <- h5_raw_path$`Gene Expression`
message(Sys.time() - st)

# Accessing raw_data with 10X Genomics experiment. Creates a SingleCellExperiment from the CellRanger
# It is +5 times faster than Read10X_h5, but it does not perform well the ranks; it needs to be loaded with Read10X_h5 to make it easy
#raw.sce <- read10xCounts(sample_path, col.names=TRUE)   

# Count the total number of cells in the gene expression matrix and store it in 'totalCells'
totalCells <- length(Cells(raw.sce))

# Compute barcode rank statistics and identify the knee and inflection points on the total count curve
bcRanks <- barcodeRanks(raw.sce, fit.bounds = c(10, 1e3))
# Barcode rank range.
colnames(bcRanks)
range((bcRanks$rank)) 
range((bcRanks$total))  
# Get knee value and add hundred points to perform a more stringent threshold. 
knee_lower <- metadata(bcRanks)$knee + 100
# Get the inflection point. Signs change.
infection <- metadata(bcRanks)$inflection

# tracking stats of read10xCounts against Read10X_h5 
message(paste(s_sample,knee_lower,infection))

# Run emptyDrops w/ knee + 100 to make it more strident 
# emptyDrops distinguish between droplets containing cells and ambient RNA in a droplet-based single-cell RNA experiment
st <- Sys.time()
message(Sys.time(), " Starting emptyDrops")
sce.out <- DropletUtils::emptyDrops(
    raw.sce,
    niters = 30000,       # number of iterations. It use the Monte Carlo p-value
    lower = knee_lower    # numeric scalar specifying the lower bound on the total UMI count
)
en <- Sys.time() - st
message(paste0(' Processing time: ', en))    

head.matrix(sce.out,n=5)
# DataFrame with 5 rows and 5 columns
# Total   LogProb    PValue   Limited       FDR
# <integer> <numeric> <numeric> <logical> <numeric>
# AAACAGCCAAACAACA-1         7        NA        NA        NA        NA
# AAACAGCCAAACATAG-1        14        NA        NA        NA        NA
# AAACAGCCAAACCCTA-1         4        NA        NA        NA        NA

# Save an object to a file
message(paste0(' Process completed. Saving data'))
s_file_name <- here('processed-data/02_empty_droplets_stats', paste0(s_sample,'_empty_droplets_object.rds'))
saveRDS(sce.out, file = s_file_name)
    
# Get significant TRUE cells based on the FDR cutoff
cells_FALSE <- 0
cells_FT <- 0
cells_TRUE <- 0
signif_TRUE <- 0
#  add arbitrary margins on a multidimensional array
cells_stat <- addmargins(table(Signif = sce.out$FDR <= FDR_cutoff,
              Limited = sce.out$Limited,
              useNA = "ifany"))
# Signif   FALSE   TRUE   <NA>    Sum
# FALSE   1548      0      0   1548
# TRUE     821  12056      0  12877
# <NA>       0      0 721642 721642
# Sum     2369  12056 721642 736067

# get specific values from a confusion mtx
cells_FALSE <- cells_stat[1,1]
cells_FT <- cells_stat[2,1]
cells_TRUE <- cells_stat[2,2]
signif_TRUE <- cells_stat[2,4]

# Calculate non Emptydroplets value
nonEmptydroplets <- (sce.out |> as.data.frame() |> filter(FDR < FDR_cutoff) |> summarise(n = n()))$n
# [1] 10359

# Calculate percentage of non Emptydroplets
per.nonemptydroplets <- ((nonEmptydroplets*100) / totalCells) #, digits = 4)
# [1] 1.46291
    
# Build a table with the stats applied and the outputs gotten
# Table with the cellranger-arc `gene expression statistics`
tab_stats <- matrix(c(s_sample, totalCells, knee_lower, infection, FDR_cutoff, 
                      cells_FALSE, cells_FT, cells_TRUE, signif_TRUE, 
                      nonEmptydroplets, per.nonemptydroplets), ncol=11, byrow=TRUE)
colnames(tab_stats) <- c('Sample_name','Total_cells','Knee_lower','Inflection','FDR_cutoff',
                         'Signif.FALSE','F/T','T/T','Signif.TRUE',
                         'Non.Emptydroplets', '%Emptydroplets')
tab_stats <- as.table(tab_stats)
        
#Outupt in csv format: 
    # Sample_name Total_cells Knee_lower Inflection FDR_cutoff Signif.FALSE  F.T
    # 1   1_HPC_KDM      736067        781        209      0.001         1548  821
    # T.T Signif.TRUE Non.Emptydroplets X.Emptydroplets
    # 1 12056       12877             12877       1.7494331

# Export the table to CSV for further analysis
message(paste0(' Process completed. Saving stats'))
s_file_name <- here('processed-data/02_empty_droplets_stats', paste0(s_sample,'_empty_droplets_stats.csv'))
write.csv(tab_stats, file=s_file_name, quote=FALSE, row.names=FALSE)

# Some code used as reference
# # create results dataframe
# # https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/03_droplet_qc.R
# 

# Droplet Elbow plot
s_file_name <- here('plots/02_empty_droplets_stats', paste0(s_sample,'_empty_droplets_knee_plot.png'))
#png(s_file_name)

define_theme <- function(size = 15) {
  theme_bw() +
    theme(text = element_text(size = size))
}

# Prepare data frame with additional FDR column
droplet_elbow_data <- as.data.frame(bcRanks) %>%
  mutate(FDR = sce.out$FDR)

# Define parameters
knee_meta <- metadata(bcRanks)$knee
knee_lower_label <- paste0("Knee est 'lower' (", knee_lower, ')')
second_knee_label <- paste0("Second Knee (", knee_meta, ')')
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
  # Apply theme
  define_theme() +
  theme(legend.position = "bottom")

# Save the png format
#dev.off()
message(paste0(' Process completed. Saving knee plot'))
ggsave(droplet_elbow_plot, filename = s_file_name) 

# # plot the elbow
# droplet_elbow_plot
# 
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


# > print("Reproducibility information:")
# 
# [1] "Reproducibility information:"
# > Sys.time()
# 
# [1] "2023-08-17 15:22:15 EDT"
# > proc.time()
# 
# user   system  elapsed 
# 1666.385   38.328 5500.404 
# > options(width = 120)
# 
# > session_info()
# 
# it64                  4.0.5      2020-08-30 [2] CRAN (R 4.3.1)
# bitops                 1.0-7      2021-04-24 [2] CRAN (R 4.3.1)
# cli                    3.6.1      2023-03-23 [2] CRAN (R 4.3.1)
# cluster                2.1.4      2022-08-22 [3] CRAN (R 4.3.1)
# codetools              0.2-19     2023-02-01 [3] CRAN (R 4.3.1)
# colorspace             2.1-0      2023-01-23 [2] CRAN (R 4.3.1)
# cowplot                1.1.1      2020-12-30 [2] CRAN (R 4.3.1)
# crayon                 1.5.2      2022-09-29 [2] CRAN (R 4.3.1)
# data.table             1.14.8     2023-02-17 [2] CRAN (R 4.3.1)
# DelayedArray           0.26.6     2023-07-02 [2] Bioconductor
# DelayedMatrixStats     1.22.1     2023-06-09 [2] Bioconductor
# deldir                 1.0-9      2023-05-17 [2] CRAN (R 4.3.1)
# digest                 0.6.33     2023-07-07 [2] CRAN (R 4.3.1)
# dotCall64              1.0-2      2022-10-03 [2] CRAN (R 4.3.1)
# dplyr                * 1.1.2      2023-04-20 [2] CRAN (R 4.3.1)
# dqrng                  0.3.0      2021-05-01 [2] CRAN (R 4.3.1)
# DropletUtils         * 1.20.0     2023-04-25 [2] Bioconductor
# edgeR                  3.42.4     2023-05-31 [2] Bioconductor
# ellipsis               0.3.2      2021-04-29 [2] CRAN (R 4.3.1)
# fansi                  1.0.4      2023-01-22 [2] CRAN (R 4.3.1)
# farver                 2.1.1      2022-07-06 [2] CRAN (R 4.3.1)
# fastDummies            1.7.3      2023-07-06 [1] CRAN (R 4.3.1)
# fastmap                1.1.1      2023-02-24 [2] CRAN (R 4.3.1)
# fitdistrplus           1.1-11     2023-04-25 [1] CRAN (R 4.3.0)
# future                 1.33.0     2023-07-01 [2] CRAN (R 4.3.1)
# future.apply           1.11.0     2023-05-21 [1] CRAN (R 4.3.1)
# generics               0.1.3      2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.36.1     2023-06-21 [2] Bioconductor
# GenomeInfoDbData       1.2.10     2023-07-20 [2] Bioconductor
# GenomicRanges        * 1.52.0     2023-04-25 [2] Bioconductor
# VP ggplot2              * 3.4.2      2023-08-14 [2] CRAN (R 4.3.1) (on disk 3.4.3)
# ggrepel                0.9.3      2023-02-03 [2] CRAN (R 4.3.1)
# ggridges               0.5.4      2022-09-26 [2] CRAN (R 4.3.1)
# globals                0.16.2     2022-11-21 [2] CRAN (R 4.3.1)
# glue                   1.6.2      2022-02-24 [2] CRAN (R 4.3.1)
# goftest                1.2-3      2021-10-07 [1] CRAN (R 4.3.0)
# gridExtra              2.3        2017-09-09 [2] CRAN (R 4.3.1)
# gtable                 0.3.3      2023-03-21 [2] CRAN (R 4.3.1)
# HDF5Array              1.28.1     2023-05-01 [2] Bioconductor
# hdf5r                  1.3.8      2023-01-21 [2] CRAN (R 4.3.1)
# here                 * 1.0.1      2020-12-13 [2] CRAN (R 4.3.1)
# htmltools              0.5.5      2023-03-23 [2] CRAN (R 4.3.1)
# htmlwidgets            1.6.2      2023-03-17 [2] CRAN (R 4.3.1)
# httpuv                 1.6.11     2023-05-11 [2] CRAN (R 4.3.1)
# httr                   1.4.7      2023-08-15 [1] CRAN (R 4.3.1)
# ica                    1.0-3      2022-07-08 [1] CRAN (R 4.3.0)
# igraph                 1.5.0      2023-06-16 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1     2023-06-22 [2] Bioconductor
# irlba                  2.3.5.1    2022-10-03 [2] CRAN (R 4.3.1)
# jsonlite               1.8.7      2023-06-29 [2] CRAN (R 4.3.1)
# KernSmooth             2.23-22    2023-07-10 [3] CRAN (R 4.3.1)
# later                  1.3.1      2023-05-02 [2] CRAN (R 4.3.1)
# lattice                0.21-8     2023-04-05 [3] CRAN (R 4.3.1)
# lazyeval               0.2.2      2019-03-15 [2] CRAN (R 4.3.1)
# leiden                 0.4.3      2022-09-10 [1] CRAN (R 4.3.0)
# lifecycle              1.0.3      2022-10-07 [2] CRAN (R 4.3.1)
# limma                  3.56.2     2023-06-04 [2] Bioconductor
# listenv                0.9.0      2022-12-16 [2] CRAN (R 4.3.1)
# lmtest                 0.9-40     2022-03-21 [2] CRAN (R 4.3.1)
# locfit                 1.5-9.8    2023-06-11 [2] CRAN (R 4.3.1)
# magrittr               2.0.3      2022-03-30 [2] CRAN (R 4.3.1)
# MASS                   7.3-60     2023-05-04 [3] CRAN (R 4.3.1)
# Matrix                 1.6-0      2023-07-08 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.2     2023-06-09 [2] Bioconductor
# matrixStats          * 1.0.0      2023-06-02 [2] CRAN (R 4.3.1)
# mime                   0.12       2021-09-28 [2] CRAN (R 4.3.1)
# miniUI                 0.1.1.1    2018-05-18 [2] CRAN (R 4.3.1)
# munsell                0.5.0      2018-06-12 [2] CRAN (R 4.3.1)
# nlme                   3.1-162    2023-01-31 [3] CRAN (R 4.3.1)
# parallelly             1.36.0     2023-05-26 [2] CRAN (R 4.3.1)
# patchwork              1.1.2      2022-08-19 [2] CRAN (R 4.3.1)
# pbapply                1.7-2      2023-06-27 [2] CRAN (R 4.3.1)
# pillar                 1.9.0      2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3      2019-09-22 [2] CRAN (R 4.3.1)
# plotly                 4.10.2     2023-06-03 [2] CRAN (R 4.3.1)
# plyr                   1.8.8      2022-11-11 [2] CRAN (R 4.3.1)
# png                    0.1-8      2022-11-29 [2] CRAN (R 4.3.1)
# polyclip               1.10-4     2022-10-20 [2] CRAN (R 4.3.1)
# progressr              0.14.0     2023-08-10 [1] CRAN (R 4.3.1)
# promises               1.2.0.1    2021-02-11 [2] CRAN (R 4.3.1)
# purrr                  1.0.1      2023-01-10 [2] CRAN (R 4.3.1)
# R.methodsS3            1.8.2      2022-06-13 [2] CRAN (R 4.3.1)
# R.oo                   1.25.0     2022-06-12 [2] CRAN (R 4.3.1)
# R.utils                2.12.2     2022-11-11 [2] CRAN (R 4.3.1)
# R6                     2.5.1      2021-08-19 [2] CRAN (R 4.3.1)
# ragg                   1.2.5      2023-01-12 [2] CRAN (R 4.3.1)
# RANN                   2.6.1      2019-01-08 [2] CRAN (R 4.3.1)
# RColorBrewer           1.1-3      2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11     2023-07-06 [2] CRAN (R 4.3.1)
# RcppAnnoy              0.0.21     2023-07-02 [2] CRAN (R 4.3.1)
# RcppHNSW               0.4.1      2022-07-18 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12  2023-03-27 [2] CRAN (R 4.3.1)
# reshape2               1.4.4      2020-04-09 [2] CRAN (R 4.3.1)
# reticulate             1.30       2023-06-09 [1] CRAN (R 4.3.1)
# rhdf5                  2.44.0     2023-04-25 [2] Bioconductor
# rhdf5filters           1.12.1     2023-04-30 [2] Bioconductor
# Rhdf5lib               1.22.0     2023-04-25 [2] Bioconductor
# rlang                  1.1.1      2023-04-28 [2] CRAN (R 4.3.1)
# ROCR                   1.0-11     2020-05-02 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3      2022-04-02 [2] CRAN (R 4.3.1)
# RSpectra               0.16-1     2022-04-24 [2] CRAN (R 4.3.1)
# Rtsne                  0.16       2022-04-17 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.4      2023-05-14 [2] Bioconductor
# S4Vectors            * 0.38.1     2023-05-02 [2] Bioconductor
# scales                 1.2.1      2022-08-20 [2] CRAN (R 4.3.1)
# scattermore            1.0        2023-05-03 [1] CRAN (R 4.3.0)
# sctransform            0.3.5      2022-09-21 [1] CRAN (R 4.3.0)
# scuttle                1.10.1     2023-05-02 [2] Bioconductor
# sessioninfo          * 1.2.2      2021-12-06 [2] CRAN (R 4.3.1)
# Seurat               * 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
# SeuratObject         * 4.9.9.9084 2023-05-08 [1] Github (mojaveazure/seurat-object@185884a)
# shiny                  1.7.4.1    2023-07-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0     2023-04-25 [2] Bioconductor
# sp                   * 2.0-0      2023-06-22 [2] CRAN (R 4.3.1)
# spam                   2.9-1      2022-08-07 [2] CRAN (R 4.3.1)
# sparseMatrixStats      1.12.2     2023-07-02 [2] Bioconductor
# spatstat.data          3.0-1      2023-03-12 [1] CRAN (R 4.3.0)
# spatstat.explore       3.2-1      2023-05-13 [1] CRAN (R 4.3.0)
# spatstat.geom          3.2-1      2023-05-09 [1] CRAN (R 4.3.0)
# spatstat.random        3.1-5      2023-05-11 [1] CRAN (R 4.3.0)
# spatstat.sparse        3.0-1      2023-03-12 [1] CRAN (R 4.3.0)
# spatstat.utils         3.0-3      2023-05-09 [1] CRAN (R 4.3.0)
# stringi                1.7.12     2023-01-11 [2] CRAN (R 4.3.1)
# stringr                1.5.0      2022-12-02 [2] CRAN (R 4.3.1)
# SummarizedExperiment * 1.30.2     2023-06-06 [2] Bioconductor
# survival               3.5-5      2023-03-12 [3] CRAN (R 4.3.1)
# systemfonts            1.0.4      2022-02-11 [2] CRAN (R 4.3.1)
# tensor                 1.5        2012-05-05 [1] CRAN (R 4.3.0)
# textshaping            0.3.6      2021-10-13 [2] CRAN (R 4.3.1)
# tibble                 3.2.1      2023-03-20 [2] CRAN (R 4.3.1)
# tidyr                  1.3.0      2023-01-24 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0      2022-10-10 [2] CRAN (R 4.3.1)
# utf8                   1.2.3      2023-01-31 [2] CRAN (R 4.3.1)
# uwot                   0.1.16     2023-06-29 [2] CRAN (R 4.3.1)
# vctrs                  0.6.3      2023-06-14 [2] CRAN (R 4.3.1)
# viridisLite            0.4.2      2023-05-02 [2] CRAN (R 4.3.1)
# withr                  2.5.0      2022-03-03 [2] CRAN (R 4.3.1)
# xtable                 1.8-4      2019-04-21 [2] CRAN (R 4.3.1)
# XVector                0.40.0     2023-04-25 [2] Bioconductor
# zlibbioc               1.46.0     2023-04-25 [2] Bioconductor
# zoo                    1.8-12     2023-04-13 [2] CRAN (R 4.3.1)
# 
# [1] /users/csoto/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
