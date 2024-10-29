########################################################################
## Read Seurat clusters to identify cell-types based on a custom marker gene list
## INPUT:
##      A Seurat batch corrected data with `min_cells`. Raw data is `GEX` multiome assays processed with `cellranger-count` pipeline 
##      A csv file with Gene-marker list.
##      A csv file with DGE from a Seurat Harmony/CCA data (not pseudo-bulked)
## 
## OUPUT:
##      1) A csv files with clusters cell-type identification
##
## Authors. CSC 
## Date. Feb 27th, 2024 / last md. August 2024
########################################################################

## load libraries
library(tidyverse)
library(dplyr)
library(data.table)
library(magrittr)
library(here)

here::here()

## read input arguments
args = commandArgs(trailingOnly=TRUE)
cellranger_pipe <- args[2]
## For testing:
# cellranger_pipe <- "CR_complementBarcodes"
# cellranger_pipe <- "CR_crossBarcodes"

## input validations
if (length(cellranger_pipe)) {
  message("CellRanger input: ", cellranger_pipe)
  if (cellranger_pipe=="CR_crossBarcodes") {
    inputDir <- here("processed-data", "03_pseudobulking", "cellranger_count")
    inputDir_cvs <- here("processed-data", "03_pseudobulking", "cellranger_count", "cvs_files_markers")
    processedDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellranger_count")
    cvsDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellranger_count", "cvs_files_markers")
  } else {
    inputDir <- here("processed-data", "03_pseudobulking", "cellranger_count_complement")
    inputDir_cvs <- here("processed-data", "03_pseudobulking", "cellranger_count_complement", "cvs_files_markers")
    processedDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellranger_count_complement")
    cvsDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellranger_count_complement", "cvs_files_markers")
  }  
  
} else {
  message("Input argument missed")
  message("CellRanger input: ", cellranger_pipe)
  stop()
}

## Check directories
if (!dir.exists(processedDir)) {dir.create(processedDir)}
if (!dir.exists(cvsDir)) {dir.create(cvsDir)}

## Contains marker lists 
source(here("code", "04_DiffExpr_Clustering_seurat", "remote_DGE_marker_gene_lists.R"))       # Call functions to read paths

get_seurat <- function(name) { sobj <- readRDS(name)}


#############################           Initials        ################################

## We are only using the norm count with Harmony
count_mtx_type <- 'norm_counts' 
Seurat_reduction <- 'Harmony' 
minCells <- 1 ### Minimum cells by cluster 

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts'}

## Build Seurat object base name and csv deg
if (Seurat_reduction=='CCA') {
  #Seurat_base_name <- paste0(Seurat_base_name, '_CCA_All')
} else {
  if (cellranger_pipe=="CR_crossBarcodes") {
    Seurat_base_name <- paste0(Seurat_base_name, '_Harmony_All')
    DGE_cvs_name <- paste0(Seurat_base_name, "markers.csv")
  } else { ##CR_complementBarcodes
    Seurat_base_name <- paste0(Seurat_base_name, '_Harmony_All_subset')
    DGE_cvs_name <- paste0(Seurat_base_name, "_markers.csv")
  }
}

## Validate Seurat object exists
if (length(list.files(inputDir, pattern = Seurat_base_name)==1)) {
  message("Processing: ", Seurat_base_name)
} else {
  message("Input seurat object missed!")
  stop()
}
# Ex. seurat.norm_counts_Harmony_All_subset.rds

message("Starting cell-type identification for ", Seurat_base_name)

## Load Seurat Integrated with cluster information
SeuratOBJ <- get_seurat(here(inputDir, paste0(Seurat_base_name, ".rds")))
## verification of the integration
print(table(SeuratOBJ$orig.ident))

## Gene markers lists
markers.custom = list()
markers.custom[["data_driven"]] <- get_Top50r_markers_genes_Hb()
markers.custom[["literature_base"]] <- get_erik_and_Hb_markers_genes()  
## sub-population list
#names(markers.custom$literature_base)
#names(markers.custom)

## set labels for the number of top DGE genes to pick up
prefix_name <- 'all_gm'                                    # prefix to save matched markers found in the clusters
n_slice <- 20  


#############################  Set the DGE list to parse  ################################

## Extract cluster data
md <- SeuratOBJ@meta.data %>% as.data.table

## Apply vertical format to unique cluster with number of UMIs, arranged by sample and cluster number
mdT <- md[, .N, by = c("orig.ident", "seurat_clusters")] %>%
    arrange(., orig.ident, seurat_clusters, .by_group = FALSE)
df_mdT <- as.data.frame(mdT)

message("Total cells in the Seurat object: ", sum(df_mdT$N))

cvs_name <- paste0(Seurat_base_name, '_cluster_info.csv')
## Save clustering information; e.g: seurat.data_counts_Harmony_cluster_info.csv
write.csv(df_mdT, here(cvsDir, cvs_name))

## extract unique clusters in ascending order
clusters <- unique(df_mdT$seurat_clusters)
clusters <- as.integer(levels(clusters)[as.integer(clusters)])

message('Identifing cell types for ', length(clusters),' Seurat clusters from ', Seurat_base_name)

## Read DGE cvs file for all clusters for the given sample
#DGE_cvs_name <- paste0(Seurat_base_name, '_',Seurat_reduction, '_Allmarkers_min', minCells, 'cells.csv')
DGE_cvs_name <- here(inputDir_cvs, DGE_cvs_name) 
seurat_clust <- as.data.frame(read.csv(DGE_cvs_name, header = TRUE))
head(seurat_clust, n=3)
# Seurat output from FindAllmarkers()
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster   gene
# 1     0  -1.836151 0.078 0.736         0     C_0  NPAS3

message('Parsing ', length(markers.custom), ' gene-markers lists on ', length(clusters) ,' clusters in ', Seurat_base_name, ' ', Seurat_reduction, ' reduction')


####### Parse the 20 Top DGE genes from GEX cluster against the marker genes list provided ####### 

names(markers.custom) #[1] "literature_base" "data_driven" 
idx_lst <- 0 

for (markers.lst in markers.custom) {
  # for testing: markers.lst <- markers.custom$literature_base
  # for testing: markers.lst <- markers.custom$data_driven
  
  ## Build df to save cell-types that match with the gene-marker-list
  all_gene_match <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Feature.ID", "Feature.Name", "Cluster.Adjusted.p.value", "cell-type", "cluster"))
  
  idx_lst <- idx_lst+1
  prefix_name <- paste0(names(markers.custom[idx_lst]), '_top', n_slice)
  print(paste("Searching markers for ", names(markers.lst)))
  
  ## Read x cluster and extract the top 10 genes
  for (clust in clusters) {
    # Testing: clust <- 2
    message("Parsing cluster ", as.character(clust))
    top_DGE_clust <- seurat_clust |> 
      dplyr::filter(cluster == clust, p_val_adj < 0.05) |> slice_head(n = n_slice)
    
    # get a vector with all cell-types with their marker genes
    gm_lst <- as.vector(as.list(markers.lst))
    i_pos <- 0      # reset gene-marker list position
    
    for ( gm in gm_lst ) {
      #for testing: gm <- gm_lst[[0]]
      i_pos <- i_pos+1                        # control cell type position
      # Match top10genes with the marker genes for the cell-type x 
      gene_match <- top_DGE_clust |> filter_all(any_vars(. %in% gm))
      # add matched genes to a dataframe
      if ( nrow(gene_match) > 0 ) {
        names(gene_match)[names(gene_match) == clust ] <- "Cluster.Adjusted.p.value" # rename cols to rbind
        gene_match['cell-type']  <- names(gm_lst[i_pos])
        gene_match['cluster']  <- clust
        all_gene_match <- rbind(all_gene_match, gene_match)
        message("Gene match for cell-type: ", names(gm_lst[i_pos]))
      } else {
        #message("None gene match for cell-type: ", names(gm_lst[i_pos]))
      }
    }
    
  }
  
  if (cellranger_pipe=="CR_crossBarcodes") {
    habenula_markers_cvs_name <- here(cvsDir, paste0(Seurat_base_name, '_cellTypes_', prefix_name, ".csv"))
  } else { # complement subset
    habenula_markers_cvs_name <- here(cvsDir, paste0(Seurat_base_name, '_cellTypes_', prefix_name, ".csv"))
  }
  print(paste("Printing results in ", habenula_markers_cvs_name))
  write.csv(all_gene_match, habenula_markers_cvs_name, row.names=FALSE)

}

head(all_gene_match, n=3)
# dim(all_gene_match)

message(' Cell type identification in clusters done!')


## slurm script reproducibility
# library("slurmjobs")
# slurmjobs::job_loop(
#   loops = list(cellranger_lst = c("CR_crossBarcodes", "CR_complementBarcodes")),
#   name = "01_Hb_celltypes_from_seurat_CRcount_v2",
#   cores = 2,
#   create_shell = TRUE
# )


library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


# > library("sessioninfo")
# > print('Reproducibility information:')
# [1] "Reproducibility information:"
# > # Last modification
#     > Sys.time()
# [1] "2024-02-27 20:41:05 EST"
# > #"2023-04-04 12:42:26 EDT"
#     > proc.time()
# user   system  elapsed 
# 52.926    7.828 2838.397 
# > options(width = 120)
# > session_info()
# CRAN (R 4.3.1)
# BPCells            0.1.0      2023-10-02 [1] Github (bnprks/BPCells@ac4376d)
# cellranger         1.1.0      2016-07-27 [2] CRAN (R 4.3.1)
# cli                3.6.1      2023-03-23 [2] CRAN (R 4.3.1)
# cluster            2.1.4      2022-08-22 [3] CRAN (R 4.3.1)
# codetools          0.2-19     2023-02-01 [3] CRAN (R 4.3.1)
# colorspace         2.1-0      2023-01-23 [2] CRAN (R 4.3.1)
# cowplot            1.1.1      2020-12-30 [2] CRAN (R 4.3.1)
# data.table       * 1.14.8     2023-02-17 [2] CRAN (R 4.3.1)
# deldir             1.0-9      2023-05-17 [2] CRAN (R 4.3.1)
# digest             0.6.33     2023-07-07 [2] CRAN (R 4.3.1)
# dotCall64          1.1-1      2023-11-28 [1] CRAN (R 4.3.1)
# dplyr            * 1.1.4      2023-11-17 [1] CRAN (R 4.3.1)
# ellipsis           0.3.2      2021-04-29 [2] CRAN (R 4.3.1)
# fansi              1.0.6      2023-12-08 [1] CRAN (R 4.3.1)
# fastDummies        1.7.3      2023-07-06 [1] CRAN (R 4.3.1)
# fastmap            1.1.1      2023-02-24 [2] CRAN (R 4.3.1)
# fitdistrplus       1.1-11     2023-04-25 [1] CRAN (R 4.3.1)
# forcats          * 1.0.0      2023-01-29 [2] CRAN (R 4.3.1)
# future             1.33.0     2023-07-01 [2] CRAN (R 4.3.1)
# future.apply       1.11.1     2023-12-21 [1] CRAN (R 4.3.1)
# generics           0.1.3      2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb       1.36.3     2023-09-07 [2] Bioconductor
# GenomeInfoDbData   1.2.10     2023-07-20 [2] Bioconductor
# GenomicRanges      1.52.0     2023-04-25 [2] Bioconductor
# ggplot2          * 3.5.0      2024-02-23 [1] CRAN (R 4.3.1)
# ggrepel            0.9.5      2024-01-10 [1] CRAN (R 4.3.1)
# ggridges           0.5.4      2022-09-26 [2] CRAN (R 4.3.1)
# globals            0.16.2     2022-11-21 [2] CRAN (R 4.3.1)
# glue               1.6.2      2022-02-24 [2] CRAN (R 4.3.1)
# goftest            1.2-3      2021-10-07 [1] CRAN (R 4.3.1)
# gridExtra          2.3        2017-09-09 [1] CRAN (R 4.3.1)
# gtable             0.3.4      2023-08-21 [2] CRAN (R 4.3.1)
# here             * 1.0.1      2020-12-13 [2] CRAN (R 4.3.1)
# hms                1.1.3      2023-03-21 [2] CRAN (R 4.3.1)
# htmltools          0.5.7      2023-11-03 [1] CRAN (R 4.3.1)
# htmlwidgets        1.6.2      2023-03-17 [2] CRAN (R 4.3.1)
# httpuv             1.6.14     2024-01-26 [1] CRAN (R 4.3.1)
# httr               1.4.7      2023-08-15 [2] CRAN (R 4.3.1)
# ica                1.0-3      2022-07-08 [1] CRAN (R 4.3.1)
# igraph             1.5.1      2023-08-10 [2] CRAN (R 4.3.1)
# IRanges            2.34.1     2023-06-22 [2] Bioconductor
# irlba              2.3.5.1    2022-10-03 [2] CRAN (R 4.3.1)
# jsonlite           1.8.7      2023-06-29 [2] CRAN (R 4.3.1)
# KernSmooth         2.23-22    2023-07-10 [3] CRAN (R 4.3.1)
# later              1.3.1      2023-05-02 [2] CRAN (R 4.3.1)
# lattice            0.21-8     2023-04-05 [3] CRAN (R 4.3.1)
# lazyeval           0.2.2      2019-03-15 [2] CRAN (R 4.3.1)
# leiden             0.4.3.1    2023-11-17 [1] CRAN (R 4.3.1)
# lifecycle          1.0.4      2023-11-07 [1] CRAN (R 4.3.1)
# listenv            0.9.0      2022-12-16 [2] CRAN (R 4.3.1)
# lmtest             0.9-40     2022-03-21 [2] CRAN (R 4.3.1)
# lubridate        * 1.9.3      2023-09-27 [1] CRAN (R 4.3.1)
# magrittr         * 2.0.3      2022-03-30 [2] CRAN (R 4.3.1)
# MASS               7.3-60     2023-05-04 [3] CRAN (R 4.3.1)
# Matrix             1.6-5      2024-01-11 [1] CRAN (R 4.3.1)
# matrixStats        1.2.0      2023-12-11 [1] CRAN (R 4.3.1)
# mime               0.12       2021-09-28 [2] CRAN (R 4.3.1)
# miniUI             0.1.1.1    2018-05-18 [2] CRAN (R 4.3.1)
# munsell            0.5.0      2018-06-12 [2] CRAN (R 4.3.1)
# nlme               3.1-163    2023-08-09 [3] CRAN (R 4.3.1)
# parallelly         1.36.0     2023-05-26 [2] CRAN (R 4.3.1)
# patchwork          1.1.3      2023-08-14 [2] CRAN (R 4.3.1)
# pbapply            1.7-2      2023-06-27 [2] CRAN (R 4.3.1)
# pillar             1.9.0      2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig          2.0.3      2019-09-22 [2] CRAN (R 4.3.1)
# plotly             4.10.4     2024-01-13 [1] CRAN (R 4.3.1)
# plyr               1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
# png                0.1-8      2022-11-29 [1] CRAN (R 4.3.1)
# polyclip           1.10-6     2023-09-27 [1] CRAN (R 4.3.1)
# progressr          0.14.0     2023-08-10 [1] CRAN (R 4.3.1)
# promises           1.2.1      2023-08-10 [2] CRAN (R 4.3.1)
# purrr            * 1.0.2      2023-08-10 [2] CRAN (R 4.3.1)
# R6                 2.5.1      2021-08-19 [2] CRAN (R 4.3.1)
# RANN               2.6.1      2019-01-08 [2] CRAN (R 4.3.1)
# RColorBrewer       1.1-3      2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp               1.0.11     2023-07-06 [2] CRAN (R 4.3.1)
# RcppAnnoy          0.0.21     2023-07-02 [2] CRAN (R 4.3.1)
# RcppHNSW           0.5.0      2023-09-19 [2] CRAN (R 4.3.1)
# RCurl              1.98-1.12  2023-03-27 [2] CRAN (R 4.3.1)
# readr            * 2.1.4      2023-02-10 [2] CRAN (R 4.3.1)
# readxl           * 1.4.3      2023-07-06 [2] CRAN (R 4.3.1)
# reshape2           1.4.4      2020-04-09 [2] CRAN (R 4.3.1)
# reticulate         1.35.0     2024-01-31 [1] CRAN (R 4.3.1)
# rlang              1.1.3      2024-01-10 [1] CRAN (R 4.3.1)
# ROCR               1.0-11     2020-05-02 [2] CRAN (R 4.3.1)
# rprojroot          2.0.4      2023-11-05 [1] CRAN (R 4.3.1)
# RSpectra           0.16-1     2022-04-24 [2] CRAN (R 4.3.1)
# Rtsne              0.16       2022-04-17 [2] CRAN (R 4.3.1)
# S4Vectors          0.38.2     2023-09-22 [1] Bioconductor
# scales             1.3.0      2023-11-28 [1] CRAN (R 4.3.1)
# scattermore        1.2        2023-06-12 [1] CRAN (R 4.3.1)
# sctransform        0.4.1      2023-10-19 [1] CRAN (R 4.3.1)
# sessioninfo      * 1.2.2      2021-12-06 [2] CRAN (R 4.3.1)
# Seurat             4.9.9.9067 2023-10-02 [1] Github (satijalab/seurat@99b9ded)
# SeuratObject       5.0.1      2023-11-17 [1] CRAN (R 4.3.1)
# shiny              1.8.0      2023-11-17 [1] CRAN (R 4.3.1)
# sp                 2.1-1      2023-10-16 [1] CRAN (R 4.3.1)
# spam               2.10-0     2023-10-23 [1] CRAN (R 4.3.1)
# spatstat.data      3.0-4      2024-01-15 [1] CRAN (R 4.3.1)
# spatstat.explore   3.2-5      2023-10-22 [1] CRAN (R 4.3.1)
# spatstat.geom      3.2-7      2023-10-20 [1] CRAN (R 4.3.1)
# spatstat.random    3.2-1      2023-10-21 [1] CRAN (R 4.3.1)
# spatstat.sparse    3.0-3      2023-10-24 [1] CRAN (R 4.3.1)
# spatstat.utils     3.0-4      2023-10-24 [1] CRAN (R 4.3.1)
# stringi            1.8.3      2023-12-11 [1] CRAN (R 4.3.1)
# stringr          * 1.5.1      2023-11-14 [1] CRAN (R 4.3.1)
# survival           3.5-7      2023-08-14 [3] CRAN (R 4.3.1)
# tensor             1.5        2012-05-05 [1] CRAN (R 4.3.1)
# tibble           * 3.2.1      2023-03-20 [2] CRAN (R 4.3.1)
# tidyr            * 1.3.0      2023-01-24 [2] CRAN (R 4.3.1)
# tidyselect         1.2.0      2022-10-10 [2] CRAN (R 4.3.1)
# tidyverse        * 2.0.0      2023-02-22 [2] CRAN (R 4.3.1)
# timechange         0.2.0      2023-01-11 [2] CRAN (R 4.3.1)
# tzdb               0.4.0      2023-05-12 [2] CRAN (R 4.3.1)
# utf8               1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# uwot               0.1.16     2023-06-29 [2] CRAN (R 4.3.1)
# vctrs              0.6.5      2023-12-01 [1] CRAN (R 4.3.1)
# viridisLite        0.4.2      2023-05-02 [2] CRAN (R 4.3.1)
# withr              3.0.0      2024-01-16 [1] CRAN (R 4.3.1)
# xtable             1.8-4      2019-04-21 [2] CRAN (R 4.3.1)
# XVector            0.40.0     2023-04-25 [2] Bioconductor
# zlibbioc           1.46.0     2023-04-25 [2] Bioconductor
# zoo                1.8-12     2023-04-13 [2] CRAN (R 4.3.1)
# 
# [1] /users/csoto/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library