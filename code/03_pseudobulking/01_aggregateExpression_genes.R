########################################################################

## Pseudo bulk with Aggregate expression (from Seurat) for RNA assay for CCA and Harmony reductions
## Authors. CSC/lcollado
## Date. March 5th, 2024
##
## Input: Seurat integrated object with samples S1 and S2 
## Output: Seurat pseudobulked and heatmap of the correlation over specific marker genes  
##
## NOTES: require ~20G mem

########################################################################

library('Seurat')                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
## Additional required packages for aggregation
library('multtest')
library('metap')
## Additional required packages for customize clusters
library('scCustomize')
library('magrittr')
library('tidyverse')

library(here)

here::here()

if (!packageVersion("Seurat")=='4.9.9.9060') {
    stop
    message('This pipeline was implemented with Seurat v5 and Signac v1.11+ ')
    message('You need the laterst Seurat v5 (‘4.9.9.9060’)')
    message('Current available repository on: https://satijalab.org/seurat/articles/install.html  ') }

# Check if processed_data directory exists, if not create it
if (!dir.exists(here("processed-data/03_pseudobulking/"))) {
    dir.create(here("processed-data/03_pseudobulking/"))
}
# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/03_pseudobulking/"))) {
    dir.create(here("plots/03_pseudobulking/"))
}


########################    Initials ########################  

## Select the count-mtx to merge (raw or normalized data)
count_mtx_type <- 'data_counts'
#count_mtx_type <- 'norm_counts' 
#Seurat_reduction <- 'CCA'
Seurat_reduction <- 'Harmony'

## load pre-existing Seurat
get_seurat <- function(name) {
    sobj <- readRDS(name)
    # verification of the integration
    print(table(sobj$orig.ident))
    # S1_Hb_KDM S2_Hb_KDM
    # 8178      9816
    print(head(sobj, n=2))
    return(sobj)
}


##### load pre-existing Seurat objects

## Compose Seurat object name
if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.combined.data_counts_PCA' } else { Seurat_base_name <- 'seurat.combined.norm_counts_PCA' }
if (Seurat_reduction=='CCA') {
  rds_name <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_CCA.rds'))
} else {
  rds_name <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_Harmony.rds'))
}
rds_name
# ~/seurat.combined.data_counts_PCA_Harmony.rds
# ~/seurat.combined.data_counts_PCA_CCA.rds

## Load Seurat object
SeuratOBJ <- get_seurat(rds_name)
# An object of class Seurat 
# 36601 features across 17994 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
# 3 layers present: data, counts, scale.data
# 4 dimensional reductions calculated: pca, umap.unintegrated, integrated.cca, umap

## verify object
table(SeuratOBJ$orig.ident)
# S1_Hb_KDM S2_Hb_KDM 
# 8178      9816 
head(colnames(SeuratOBJ))
tail(colnames(SeuratOBJ))
colnames(SeuratOBJ@meta.data)



###################### Pseudo bulk expression data  ######################

## Run in an integrated Seurat
SeuratOBJ[["RNA"]] <- JoinLayers(SeuratOBJ[["RNA"]])

## To avoid issue when having few cells need to adjust the minimum number of cells
## For example, if there is a cluster "15" that has 0 cells, the function will skip that cluster with a warning (that's perfect). Also if the number of cells is between min.cells.groups (default = 1-3) and 0, an error is thrown and it stops working. That is why I previously remove from the Seurat Object the cells of the clusters with 3 or less cells for each condition/sample. 
few_cells_samples <- unique(SeuratOBJ@meta.data$orig.ident)
few_cells <- vector()

# Calculate the minimum number of cells by cluster. In this case is set to 1
# library(arsenal)
# sclust <- as.data.frame(SeuratOBJ@meta.data$seurat_clusters)
# sclust_r <- as.data.frame(SeuratOBJ@meta.data$seurat_clusters.renamed)
# comparedf(sclust, sclust_r)

for (i in 1:length(few_cells_samples)) {   # remove cellstype w/ less than 1 cells in each sample/condition
  # for testing: i <- 1
  few_cells_tmp <- table(SeuratOBJ@meta.data$seurat_clusters[SeuratOBJ@meta.data$orig.ident == few_cells_samples[i]]) <= 1
  few_cells_tmp <- names(few_cells_tmp)[few_cells_tmp == "TRUE"]
  few_cells <- c(few_cells,few_cells_tmp)
}

message(' Number of clusters with >1 cell: ', length(few_cells), '; Cluster ID:', few_cells)
# > few_cells
# [1] "18"

clusters <- sort(unique(SeuratOBJ@meta.data$seurat_clusters))
`%notin%` <- Negate(`%in%`) 
clusters <- clusters[clusters %notin% few_cells]

# count cells by clusters
SeuratOBJ@meta.data$seurat_clusters
table(Idents(SeuratOBJ))

## Set target clusters and genes to aggregate
## Note that we previously search Habenula marker genes at the top20 DEG for each cluster (out ~/.processed-data/04_DiffExpr_Clustering_seurat/)
##      - Marker genes comes from Hb_pilot project 

if (Seurat_reduction=='CCA') {
  clust_selected <- c(8, 15, 18)
  ## genes selected from Hb_pilot, defined as Hb general, LHb and MHb
  markers.to.plot <- c("MMRN1", "GPR151", "POU4F1", # Hb
                       "LINC01876", # LHb
                       "CD24", "AC004594.1") # MHb
} else {
  clust_selected <- c(6, 12, 14, 15)
  ## genes selected from Hb_pilot, defined as Hb general, LHb and MHb
  markers.to.plot <- c("MMRN1", "GPR151", "POU4F1", # Hb
                       "EPHA5", "TLL1", # LHb
                       "AC109466.1", "AC008415.1", "GPR149", "GNG8", "NEUROD1", 
                       "RASGRP1", "SLC5A7", "CHRNB3", "SCUBE1", "LINC02143", "CD24") # MHb  
}  


## Subset the clusters and marker genes selected
SeuratOBJ_Hb_all <- subset(SeuratOBJ, 
                           subset = seurat_clusters %in% clust_selected)

## look at minimum cells by cluster assigned
clust_tab <- table(SeuratOBJ_Hb_all$seurat_clusters)
min(clust_tab[clust_tab > 0])
# S1:
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
# 0   0   0   0   0   0 830   0   0   0   0   0 180   0 102  99   0 

## look at minimum cells by cluster assigned by sample
md <- SeuratOBJ_Hb_all@meta.data %>% as.data.table
## Apply vertical format to unique cluster with number of umis, arranged by sample and cluster number
mdT <- md[, .N, by = c("orig.ident", "seurat_clusters")] %>%
  arrange(., orig.ident, seurat_clusters, .by_group = FALSE)
df_mdT <- as.data.frame(mdT)
df_mdT 
# S1: orig.ident seurat_clusters   N
# 1  S1_Hb_KDM               6 665
# 2  S1_Hb_KDM              12 133
# 3  S1_Hb_KDM              14  29
# 4  S1_Hb_KDM              15  13
# 5  S2_Hb_KDM               6 165
# 6  S2_Hb_KDM              12  47
# 7  S2_Hb_KDM              14  73
# 8  S2_Hb_KDM              15  86


## perform pseudo bulk 
# If return.seurat = TRUE, aggregated values are placed in the 'counts' layer of the returned object. The data is then normalized by running NormalizeData on the aggregated counts. ScaleData is then run on the default assay before returning the object.
# Defaults to: normalization.method = "LogNormalize", scale.factor = 10000

SeuratOBJ_Hb_all_pseudobulked <- AggregateExpression(SeuratOBJ_Hb_all, return.seurat = TRUE, 
                                                     group.by = c("seurat_clusters", "orig.ident"))
# S1: An object of class Seurat 
# 36601 features across 8 samples within 1 assay 
# Active assay: RNA (36601 features, 0 variable features)
# 3 layers present: counts, data, scale.data
table(SeuratOBJ_Hb_all_pseudobulked$orig.ident)
table(SeuratOBJ_Hb_all_pseudobulked$seurat_clusters)


## plot the pseudo bulk 

pdf_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_DoHeatmap_pseudobulk.pdf')
pdf_name <- here('plots/04_DiffExpr_Clustering_seurat', pdf_file)
pdf(file = pdf_name)
# ~/seurat.combined.data_counts_PCA_Harmony_DoHeatmap_pseudobulk.pdf

p1 <- DoHeatmap(object = SeuratOBJ_Hb_all_pseudobulked, 
                features=markers.to.plot, label = TRUE, angle=45, group.by = "seurat_clusters",
                size=4) 
#p1 + scale_fill_gradientn(limits = c(-2, 2), colours = PurpleAndYellow(), na.value = "white")

dev.off()


## Save new Seurat pseudo bulk 

rds_name <- paste0(Seurat_base_name,'_', Seurat_reduction, '_pseudobulk.rds')
rds_name <- here('processed-data/03_pseudobulking', rds_name)
# file name: ~/.seurat.combined.data_counts_PCA_Harmony_pseudobulk.rds
saveRDS(SeuratOBJ, file = rds_name)
message('Seurat pseudobulk completed! ')   



############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()


# > library("sessioninfo")
# > print('Reproducibility information:')
# [1] "Reproducibility information:"
# > # Last modification
#   > Sys.time()
# [1] "2024-03-13 13:46:47 EDT"
# > proc.time()
# user   system  elapsed 
# 354.683   20.761 1747.221 
# > options(width = 120)
# > session_info()
# 
# cowplot            1.1.3      2024-01-22 [2] CRAN (R 4.3.2)
# data.table       * 1.15.0     2024-01-30 [2] CRAN (R 4.3.2)
# deldir             2.0-2      2023-11-23 [2] CRAN (R 4.3.2)
# digest             0.6.34     2024-01-11 [2] CRAN (R 4.3.2)
# dotCall64          1.1-1      2023-11-28 [2] CRAN (R 4.3.2)
# dplyr            * 1.1.4      2023-11-17 [2] CRAN (R 4.3.2)
# ellipsis           0.3.2      2021-04-29 [2] CRAN (R 4.3.2)
# fansi              1.0.6      2023-12-08 [2] CRAN (R 4.3.2)
# farver             2.1.1      2022-07-06 [2] CRAN (R 4.3.2)
# fastDummies        1.7.3      2023-07-06 [2] CRAN (R 4.3.2)
# fastmap            1.1.1      2023-02-24 [2] CRAN (R 4.3.2)
# fitdistrplus       1.1-11     2023-04-25 [2] CRAN (R 4.3.2)
# forcats          * 1.0.0      2023-01-29 [2] CRAN (R 4.3.2)
# future             1.33.1     2023-12-22 [2] CRAN (R 4.3.2)
# future.apply       1.11.1     2023-12-21 [2] CRAN (R 4.3.2)
# generics           0.1.3      2022-07-05 [2] CRAN (R 4.3.2)
# ggbeeswarm         0.7.2      2023-04-29 [2] CRAN (R 4.3.2)
# ggplot2          * 3.4.4      2023-10-12 [2] CRAN (R 4.3.2)
# ggprism            1.0.4      2022-11-04 [1] CRAN (R 4.3.2)
# ggrastr            1.0.2      2023-06-01 [2] CRAN (R 4.3.2)
# ggrepel            0.9.5      2024-01-10 [2] CRAN (R 4.3.2)
# ggridges           0.5.6      2024-01-23 [2] CRAN (R 4.3.2)
# GlobalOptions      0.1.2      2020-06-10 [2] CRAN (R 4.3.2)
# globals            0.16.2     2022-11-21 [2] CRAN (R 4.3.2)
# glue               1.7.0      2024-01-09 [2] CRAN (R 4.3.2)
# goftest            1.2-3      2021-10-07 [2] CRAN (R 4.3.2)
# gridExtra          2.3        2017-09-09 [2] CRAN (R 4.3.2)
# gtable             0.3.4      2023-08-21 [2] CRAN (R 4.3.2)
# here             * 1.0.1      2020-12-13 [2] CRAN (R 4.3.2)
# hms                1.1.3      2023-03-21 [2] CRAN (R 4.3.2)
# htmltools          0.5.7      2023-11-03 [2] CRAN (R 4.3.2)
# htmlwidgets        1.6.4      2023-12-06 [2] CRAN (R 4.3.2)
# httpuv             1.6.14     2024-01-26 [2] CRAN (R 4.3.2)
# httr               1.4.7      2023-08-15 [2] CRAN (R 4.3.2)
# ica                1.0-3      2022-07-08 [2] CRAN (R 4.3.2)
# igraph             2.0.1.9008 2024-02-09 [2] Github (igraph/rigraph@39158c6)
# irlba              2.3.5.1    2022-10-03 [2] CRAN (R 4.3.2)
# janitor            2.2.0      2023-02-02 [1] CRAN (R 4.3.2)
# jsonlite           1.8.8      2023-12-04 [2] CRAN (R 4.3.2)
# KernSmooth         2.23-22    2023-07-10 [3] CRAN (R 4.3.2)
# labeling           0.4.3      2023-08-29 [2] CRAN (R 4.3.2)
# later              1.3.2      2023-12-06 [2] CRAN (R 4.3.2)
# lattice            0.22-5     2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval           0.2.2      2019-03-15 [2] CRAN (R 4.3.2)
# leiden             0.4.3.1    2023-11-17 [2] CRAN (R 4.3.2)
# lifecycle          1.0.4      2023-11-07 [2] CRAN (R 4.3.2)
# listenv            0.9.1      2024-01-29 [2] CRAN (R 4.3.2)
# lmtest             0.9-40     2022-03-21 [2] CRAN (R 4.3.2)
# lubridate        * 1.9.3      2023-09-27 [2] CRAN (R 4.3.2)
# magrittr         * 2.0.3      2022-03-30 [2] CRAN (R 4.3.2)
# MASS               7.3-60.0.1 2024-01-13 [3] CRAN (R 4.3.2)
# mathjaxr           1.6-0      2022-02-28 [1] CRAN (R 4.3.2)
# Matrix             1.6-5      2024-01-11 [3] CRAN (R 4.3.2)
# matrixStats        1.2.0      2023-12-11 [2] CRAN (R 4.3.2)
# metap            * 1.9        2023-10-09 [1] CRAN (R 4.3.2)
# mime               0.12       2021-09-28 [2] CRAN (R 4.3.2)
# miniUI             0.1.1.1    2018-05-18 [2] CRAN (R 4.3.2)
# mnormt             2.1.1      2022-09-26 [1] CRAN (R 4.3.2)
# multcomp           1.4-25     2023-06-20 [2] CRAN (R 4.3.2)
# multtest         * 2.58.0     2023-10-24 [2] Bioconductor
# munsell            0.5.0      2018-06-12 [2] CRAN (R 4.3.2)
# mutoss             0.1-13     2023-03-14 [1] CRAN (R 4.3.2)
# mvtnorm            1.2-4      2023-11-27 [2] CRAN (R 4.3.2)
# nlme               3.1-164    2023-11-27 [3] CRAN (R 4.3.2)
# numDeriv           2016.8-1.1 2019-06-06 [2] CRAN (R 4.3.2)
# paletteer          1.6.0      2024-01-21 [2] CRAN (R 4.3.2)
# parallelly         1.36.0     2023-05-26 [2] CRAN (R 4.3.2)
# patchwork          1.2.0      2024-01-08 [2] CRAN (R 4.3.2)
# pbapply            1.7-2      2023-06-27 [2] CRAN (R 4.3.2)
# pillar             1.9.0      2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig          2.0.3      2019-09-22 [2] CRAN (R 4.3.2)
# plotly             4.10.4     2024-01-13 [2] CRAN (R 4.3.2)
# plotrix            3.8-4      2023-11-10 [2] CRAN (R 4.3.2)
# plyr               1.8.9      2023-10-02 [2] CRAN (R 4.3.2)
# png                0.1-8      2022-11-29 [2] CRAN (R 4.3.2)
# polyclip           1.10-6     2023-09-27 [2] CRAN (R 4.3.2)
# presto           * 1.0.0      2024-03-08 [1] Github (immunogenomics/presto@31dc97f)
# progressr          0.14.0     2023-08-10 [2] CRAN (R 4.3.2)
# promises           1.2.1      2023-08-10 [2] CRAN (R 4.3.2)
# purrr            * 1.0.2      2023-08-10 [2] CRAN (R 4.3.2)
# qqconf             1.3.2      2023-04-14 [1] CRAN (R 4.3.2)
# R6                 2.5.1      2021-08-19 [2] CRAN (R 4.3.2)
# RANN               2.6.1      2019-01-08 [2] CRAN (R 4.3.2)
# rbibutils          2.2.16     2023-10-25 [2] CRAN (R 4.3.2)
# RColorBrewer       1.1-3      2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp             * 1.0.12     2024-01-09 [2] CRAN (R 4.3.2)
# RcppAnnoy          0.0.22     2024-01-23 [2] CRAN (R 4.3.2)
# RcppHNSW           0.6.0      2024-02-04 [2] CRAN (R 4.3.2)
# Rdpack             2.6        2023-11-08 [2] CRAN (R 4.3.2)
# readr            * 2.1.5      2024-01-10 [2] CRAN (R 4.3.2)
# rematch2           2.1.2      2020-05-01 [2] CRAN (R 4.3.2)
# reshape2           1.4.4      2020-04-09 [2] CRAN (R 4.3.2)
# reticulate         1.35.0     2024-01-31 [2] CRAN (R 4.3.2)
# rlang              1.1.3      2024-01-10 [2] CRAN (R 4.3.2)
# ROCR               1.0-11     2020-05-02 [2] CRAN (R 4.3.2)
# rprojroot          2.0.4      2023-11-05 [2] CRAN (R 4.3.2)
# RSpectra           0.16-1     2022-04-24 [2] CRAN (R 4.3.2)
# Rtsne              0.17       2023-12-07 [2] CRAN (R 4.3.2)
# sandwich           3.1-0      2023-12-11 [2] CRAN (R 4.3.2)
# scales             1.3.0      2023-11-28 [2] CRAN (R 4.3.2)
# scattermore        1.2        2023-06-12 [2] CRAN (R 4.3.2)
# scCustomize      * 2.1.2      2024-02-28 [1] CRAN (R 4.3.2)
# sctransform        0.4.1      2023-10-19 [2] CRAN (R 4.3.2)
# sessioninfo      * 1.2.2      2021-12-06 [2] CRAN (R 4.3.2)
# Seurat           * 5.0.1      2023-11-17 [2] CRAN (R 4.3.2)
# SeuratObject     * 5.0.1      2023-11-17 [2] CRAN (R 4.3.2)
# shape              1.4.6      2021-05-19 [2] CRAN (R 4.3.2)
# shiny              1.8.0      2023-11-17 [2] CRAN (R 4.3.2)
# sn                 2.1.1      2023-04-04 [1] CRAN (R 4.3.2)
# snakecase          0.11.1     2023-08-27 [1] CRAN (R 4.3.2)
# sp               * 2.1-3      2024-01-30 [2] CRAN (R 4.3.2)
# spam               2.10-0     2023-10-23 [2] CRAN (R 4.3.2)
# spatstat.data      3.0-4      2024-01-15 [2] CRAN (R 4.3.2)
# spatstat.explore   3.2-6      2024-02-01 [2] CRAN (R 4.3.2)
# spatstat.geom      3.2-8      2024-01-26 [2] CRAN (R 4.3.2)
# spatstat.random    3.2-2      2023-11-29 [2] CRAN (R 4.3.2)
# spatstat.sparse    3.0-3      2023-10-24 [2] CRAN (R 4.3.2)
# spatstat.utils     3.0-4      2023-10-24 [2] CRAN (R 4.3.2)
# stringi            1.8.3      2023-12-11 [2] CRAN (R 4.3.2)
# stringr          * 1.5.1      2023-11-14 [2] CRAN (R 4.3.2)
# survival           3.5-7      2023-08-14 [3] CRAN (R 4.3.2)
# tensor             1.5        2012-05-05 [2] CRAN (R 4.3.2)
# TFisher            0.2.0      2018-03-21 [1] CRAN (R 4.3.2)
# TH.data            1.1-2      2023-04-17 [2] CRAN (R 4.3.2)
# tibble           * 3.2.1      2023-03-20 [2] CRAN (R 4.3.2)
# tidyr            * 1.3.1      2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect         1.2.0      2022-10-10 [2] CRAN (R 4.3.2)
# tidyverse        * 2.0.0      2023-02-22 [2] CRAN (R 4.3.2)
# timechange         0.3.0      2024-01-18 [2] CRAN (R 4.3.2)
# tzdb               0.4.0      2023-05-12 [2] CRAN (R 4.3.2)
# utf8               1.2.4      2023-10-22 [2] CRAN (R 4.3.2)
# uwot               0.1.16     2023-06-29 [2] CRAN (R 4.3.2)
# vctrs              0.6.5      2023-12-01 [2] CRAN (R 4.3.2)
# vipor              0.4.7      2023-12-18 [2] CRAN (R 4.3.2)
# viridisLite        0.4.2      2023-05-02 [2] CRAN (R 4.3.2)
# withr              3.0.0      2024-01-16 [2] CRAN (R 4.3.2)
# xtable             1.8-4      2019-04-21 [2] CRAN (R 4.3.2)
# zoo                1.8-12     2023-04-13 [2] CRAN (R 4.3.2)
# 
# [1] /users/csoto/R/4.3.x
# [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library