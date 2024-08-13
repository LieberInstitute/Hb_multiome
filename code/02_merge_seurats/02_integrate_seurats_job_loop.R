########################################################################
## Merge Seurat objects to prepare for batch effect correction
##
## Authors. CSC
## Date. Feb 22, 2024
## Last modification. Aug 2024
##
## Input: Seurat RDS merged for integration 
## Output:  Seurat object with new reduction slot: CCA / Harmony
##          Plots for reference after correction    
##
## NOTES: 30G free mem recommended for 20K cells
## For slurm env: runsrun --x11 --pty --partition=interactive bash
########################################################################

## slurm script reproducibility

# slurmjobs::job_loop(
#   loops = list(type_mtx = c("data_counts", "norm_counts")),
#   name = "02_integrate_seurats_job_loop.R",
#   cores = 2,
#   create_shell = TRUE
# )

library(Seurat)                                
#library(Signac)                                
library(harmony)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(here)

here::here()

## Directory to save variable features 
dir_csv <- here("processed-data", "02_merge_seurats", "csv_files")
if (!dir.exists(dir_csv)) dir.create(dir_csv)
plotsDir <- here("plots", "02_merge_seurats")


########################    Initials ########################  

## Available for Seurat data_counts and norm_counts
##        count_mtx_type_label <- 'data_counts'      
##        count_mtx_type_label <- 'norm_counts' 

args = commandArgs(trailingOnly=TRUE)
## read count mtx type (abs_counts and normalized_counts)
count_mtx_type_label <- args[2]
## For testing:
# count_mtx_type_label <- 'data_counts'      
# count_mtx_type_label <- 'norm_counts' 

message('\nIntegrating Seurat objects for `', count_mtx_type_label, '` assays')

## RDS suffix name and input assay to process 
if (count_mtx_type_label=='data_counts') { 
  #s_sample <- 'seurat.integrated.data_counts' 
  s_sample <- 'seurat.unintegrated.data_counts' 
  s_pattern <- 'seurat.combined.data_counts' 
  rna_layer <- 'counts'
} else if (count_mtx_type_label=='norm_counts') { # normalized counts
  #s_sample <- 'seurat.integrated.norm_counts' 
  s_sample <- 'seurat.unintegrated.norm_counts' 
  s_pattern <- 'seurat.combined.norm_counts' 
  rna_layer <- 'data'
} else {
  stop()
}

message('\nSuffix name used to save unintegrated seurat is `', s_sample, '`')
## Note unintegrated data are only combined 

## read directory with Seurat objects
processedDir <- here('processed-data', '02_merge_seurats')
all_rds <- paste0(processedDir, '/', list.files(processedDir, pattern= paste0(s_pattern,".rds"))) #, recursive = TRUE
all_rds

## plot reductions calculated: pca, umpa, CCA and Harmony
plot_clust <- function(sobj, f_name, reduct, ga2) {
    
    # integrate the samples and clusters
    p1 <- DimPlot(sobj, 
                  reduction = reduct, group.by = c("orig.ident", ga2))
    png_name <- here(plotsDir, paste0(f_name,'_dimplot.png'))  
    ggsave(p1, filename = png_name, height = 5, width = 10)
    
    # visualize the two conditions side-by-side
    p1 <- DimPlot(sobj, 
                  reduction = reduct, split.by = "orig.ident")
    png_name <- here(plotsDir, paste0(f_name,'_dimplot_splitted.png'))  
    ggsave(p1, filename = png_name, height = 5, width = 10)
    
}

## Save variable features before correction at 10, 20, 50, 100 and all
save_VFeatures <- function(sobj, f_name, suffix) {
    
    # most highly variable genes
    VF <- c(10,20,50,100, 500, 1000)
    
    for (x in VF) {
        top <- head(VariableFeatures(sobj), x)
        cvs_name <- here(processedDir, "csv_files", paste0(suffix, '_', f_name, "_", as.character(x), 'VariableFeatures.csv'))
        print(cvs_name)
        write.csv(top, file.path(cvs_name), row.names=FALSE)
    }

}



############  Process PCA

message('Processing PCA for `', s_sample, '`\n')

SeuratOBJ.combined <-readRDS(all_rds[1])
## Exploration
table(SeuratOBJ.combined$orig.ident)
# 4S_Hb_KDM 5S_Hb_KDM 6S_Hb_KDM 
# 4513      1720      4335
# head(colnames(SeuratOBJ.combined))
# tail(colnames(SeuratOBJ.combined))
# colnames(SeuratOBJ.combined@meta.data)



######### Perform analysis in unintegrated data

## split the RNA measurements into two layers one for each sample

tryCatch( {
  SeuratOBJ.combined[["RNA"]] <- split(SeuratOBJ.combined[["RNA"]], f = SeuratOBJ.combined$orig.ident)
return(s) }
  , error = function(e) {print('layers are already split') } )
# Note Warning: Assay RNA changing from Assay to Assay5


## Run standard analysis: Calculate PCA cell embeddings
SeuratOBJ.combined <- NormalizeData(SeuratOBJ.combined)
## Identifies outliers on a 'mean variability plot'.
## vst method fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess).
SeuratOBJ.combined <- FindVariableFeatures(SeuratOBJ.combined, selection.method = "vst") 
## Save the `Variable Features` at top 10,20,50 and 100 
save_VFeatures(SeuratOBJ.combined, 'pca', s_sample)
## Scale data. Perform “LogNormalize” method to the GEX for each cell by the total expression multiply by a scale factor (10,000 by default), and log-transforms the result
all.genes <- rownames(SeuratOBJ.combined)
SeuratOBJ.combined <- ScaleData(SeuratOBJ.combined, features = all.genes)
## Run PCA
SeuratOBJ.combined <- RunPCA(SeuratOBJ.combined)
#Reductions(SeuratOBJ.combined)
#SeuratOBJ.combined@reductions

p1 <- ElbowPlot(SeuratOBJ.combined)
ggsave(p1, filename = here(plotsDir, paste0(s_sample, '_pca_elbow.png')), height = 4, width = 5)

## Compute nearest neighbor graph + SNN and Lovain for find clusters
SeuratOBJ.combined <- FindNeighbors(SeuratOBJ.combined, dims = 1:30, reduction = "pca")
SeuratOBJ.combined <- FindClusters(SeuratOBJ.combined, 
                          resolution = 1, cluster.name = "unintegrated_clusters")
#unique(SeuratOBJ.combined$seurat_clusters)
#unique(SeuratOBJ.combined$unintegrated_clusters)

SeuratOBJ.combined <- RunUMAP(SeuratOBJ.combined, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
head(SeuratOBJ.combined, n=2)
#SeuratOBJ.combined@reductions

# create and save UMAP-PCA plots grouped by sample and clusters and split side-by-side
#plot_clust(SeuratOBJ.combined, paste0(s_sample, '_umap.unintegrated_'), 'umap.unintegrated', 'seurat_clusters')
plot_clust(SeuratOBJ.combined, paste0(s_sample, '_umap'), 'umap.unintegrated', 'seurat_clusters')

# visualize more variable features in a heatmap
p1 <- DimHeatmap(SeuratOBJ.combined, reduction = 'pca', nfeatures = 30, fast = FALSE)
#ggsave(p1, filename = here(plotsDir, paste0(s_sample, '_umap.unintegrated_pca_heatmap.png'))) 
ggsave(p1, filename = here(plotsDir, paste0(s_sample, '_pca_heatmap.png')), height = 5, width = 5) 

# Save RDS Object
rds_name <- here(processedDir, paste0(s_sample, '_PCA.rds'))
# .../seurat.combined.data_counts_PCA.rds
saveRDS(SeuratOBJ.combined, file = rds_name)
message('Seurat unintegrated saved in ', rds_name)   



###################################################################### 
#####       CCAIntegration Integration                           ##### 
###################################################################### 
## Anchor-based CCA integration (method=CCAIntegration)
## CCAIntegration: https://satijalab.org/seurat/reference/ccaintegration
## “Seurat CCA” has the assumption that biologically more similar cells from different batches have a higher mathematical similarity (i.e. the dot product), and similarly, MNN assume similar cells from different batches have smaller Euclidean distance defined in the algorithm.

## Rename samples 
if (count_mtx_type_label=='data_counts') { 
  #s_sample <- 'seurat.integrated.data_counts'
  s_sample <- 'seurat.data_counts'
} else if (count_mtx_type_label=='norm_counts') { # normalized counts
  #s_sample <- 'seurat.integrated.norm_counts' 
  s_sample <- 'seurat.norm_counts' 
}

### Start from here IF you have pre-existing combined Seurat  with PCA
# if (count_mtx_type=='data_counts') { s_sample <- 'seurat.combined.data_counts' 
# } else { s_sample <- 'seurat.combined.norm_counts' }
# rds_name <- here(processedDir, paste0(s_sample, '_PCA.rds'))   # Seurat.combined.raw_PCA.rds
# SeuratOBJ.combined <- get_seurat(rds_name)

table(SeuratOBJ.combined$`orig.ident`)
# 4S_Hb_KDM 5S_Hb_KDM 6S_Hb_KDM 
# 4513      1720      4335

message("\nRunning Seurat-CCA Integration - ", Sys.time())

integration_method = 'CCAIntegration'
reduction_name <- 'integrated.cca'

seed(13082024)
SeuratOBJ.combined <- IntegrateLayers(object = SeuratOBJ.combined, # Default dims: 1:30
                             method = CCAIntegration, 
                             orig.reduction = "pca", 
                             new.reduction = reduction_name,
                             verbose = FALSE)

# Note. specify `k.anchor` to increase the strength of integration, add:  k.anchor = 20
# At point you have:
#   layers present: counts, data and scale.data
#   dimensional reductions calculated: pca, umap.unintegrated, integrated.cca
colnames(head(SeuratOBJ.combined))

## re-join layers after integration
SeuratOBJ.combined[["RNA"]] <- JoinLayers(SeuratOBJ.combined[["RNA"]])

message("Seurat-CCA Integration Done! - ", Sys.time())


## Cluster based on `integrated.cca` dimensionality reduction

#save_VFeatures(SeuratOBJ.combined, 'integrated.cca')
SeuratOBJ.combined <- FindNeighbors(SeuratOBJ.combined, reduction = reduction_name, dims = 1:30)
SeuratOBJ.combined <- FindClusters(SeuratOBJ.combined, resolution = 1)  # shared nearest neighbor (SNN)
# Note that 'seurat_clusters' will be overwritten everytime FindClusters is run

#SeuratOBJ.combined@reductions$integrated.cca
# A dimensional reduction object with key integratedcca_ 
table(SeuratOBJ.combined$orig.ident)
table(Idents(SeuratOBJ.combined))
#head(SeuratOBJ.combined, n=2)

## Run umap in the corresponding reduction
SeuratOBJ.combined <- RunUMAP(SeuratOBJ.combined, dims = 1:30, reduction = reduction_name)

## create and save UMAP-integrated.cca/harmony plots grouped by sample and clusters split side-by-side
plot_clust(SeuratOBJ.combined, paste0(s_sample, '_', reduction_name, '_umap'), 'umap', 'seurat_clusters')

## visualize most variable features in a heatmap
p1 <- DimHeatmap(SeuratOBJ.combined, reduction = reduction_name, nfeatures = 30, fast = FALSE)
png_name <- 
ggsave(p1, filename = here(plotsDir, paste0(s_sample, '_', reduction_name, '_heatmap.png')))#, height = 10, width = 10

# saveRDS(SeuratOBJ.combined, file = here(processedDir, paste0(s_sample, '_PCA_CCA.rds')))
saveRDS(SeuratOBJ.combined, file = here(processedDir, paste0(s_sample, '_CCA.rds')))
message('Seurat combined saved in ', rds_name)   




###################################################################### 
#####       Harmony Integration method                          ##### 
###################################################################### 
# HarmonyIntegration: https://satijalab.org/seurat/reference/harmonyintegration

## NOTE: if correction is handled in the same Seurat integrated, then Seurat Clusters are overwriting, 
##       that is why the output is saved as new object

### Start from here / load pre-existing Seurat combined objects
# if (count_mtx_type=='data_counts') { s_sample <- 'seurat.combined.data_counts' } else { s_sample <- 'seurat.combined.norm_counts' }
# rds_name <- here('processed-data/02_merge_seurats', paste0(s_sample, '_PCA.rds'))   # eurat.combined.raw_PCA.rds
# SeuratOBJ.combined <- get_seurat(rds_name)
# table(SeuratOBJ.combined$`orig.ident`)
# split the RNA measurements into two layers one for each sample, because you need to have a v5 assay (a bug ??)
#SeuratOBJ.combined[["RNA"]] <- split(SeuratOBJ.combined[["RNA"]], f = SeuratOBJ.combined$orig.ident)   #It is converted to v5 assay

message("Running Seurat-Harmony Integration - ", Sys.time())

integration_method = 'harmony'
reduction_name <- 'integrated.harmony'

seed(13082024)
# run correction. 
# max_iter=10 and up to 10 correction steps are expected. However, early_stop=TRUE so harmony will stop after the cost plateaus.
# Returns an object with a new dimensionality reduction
SeuratOBJ.combined <- SeuratOBJ.combined %>%
    RunHarmony(group.by.vars = "orig.ident",
               reduction = "pca",
               assay.use = 'RNA',
               reduction.save = reduction_name,
               plot_convergence = FALSE,
               nclust = 50,                     # Number of clusters in model. nclust=1 equivalent to simple linear regression
               max.iter = 10,                   # One round of Harmony involves one clustering and one correction step
               #max.iter.cluster = 20,          # Maximum number of rounds to run clustering at each round of Harmony
               early_stop = T
               #dims.use = 30
               )
# Harmony converged after `x` iterations
#SeuratOBJ.combined@reductions

## re-join layers after integration
SeuratOBJ.combined[["RNA"]] <- JoinLayers(SeuratOBJ.combined[["RNA"]])

message("Finishing Seurat-Harmony Integration - ", Sys.time())

## Cluster based in the new reduction
SeuratOBJ.combined <- FindNeighbors(SeuratOBJ.combined, reduction = reduction_name, dims = 1:30)
SeuratOBJ.combined <- FindClusters(SeuratOBJ.combined, resolution = 1)
# 1 singletons identified. 17 final clusters.
#SeuratOBJ.combined@reductions$integrated.harmony
table(SeuratOBJ.combined$orig.ident)
table(Idents(SeuratOBJ.combined))
#head(SeuratOBJ.combined, n=2)

# Visualization
SeuratOBJ.combined <- RunUMAP(SeuratOBJ.combined, dims = 1:30, reduction = reduction_name)

# create and save UMAP-integrated.cca plots grouped by sample and clusters and splitted side-by-side
#plot_clust(SeuratOBJ.combined, paste0(s_sample, '_umap.', reduction_name), 'umap', 'seurat_clusters')
plot_clust(SeuratOBJ.combined, paste0(s_sample, '_', reduction_name, '_umap'), 'umap', 'seurat_clusters')

# visualize more variable features in a heatmap
p1 <- DimHeatmap(SeuratOBJ.combined, reduction = reduction_name, nfeatures = 30, fast = FALSE)
png_name <- here(plotsDir, paste0(s_sample, '_', reduction_name, '_heatmap.png'))
ggsave(p1, filename = png_name)

#rds_name <- here(processedDir, paste0(s_sample, '_PCA_Harmony.rds'))
saveRDS(SeuratOBJ.combined, file = here(processedDir, paste0(s_sample, '_Harmony.rds')))
message('Seurat combined saved in ', rds_name)   


# INTEGRATION methods for Seurat V5:  https://satijalab.org/seurat/articles/seurat5_integration (Oct 31, 2023)
# https://satijalab.org/seurat/articles/integration_introduction.html (Nov 16, 2023)





############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()


# > Sys.time()
# [1] "2024-03-07 14:21:12 EST"
# > #"2023-04-04 12:42:26 EDT"
#   > proc.time()
# user   system  elapsed 
# 1065.858   74.277 3232.752 
# > options(width = 120)
# > session_info()
# 39m CRAN (R 4.3.2)
# data.table         1.15.0     2024-01-30 [2] CRAN (R 4.3.2)
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
# ggplot2          * 3.4.4      2023-10-12 [2] CRAN (R 4.3.2)
# ggrepel            0.9.5      2024-01-10 [2] CRAN (R 4.3.2)
# ggridges           0.5.6      2024-01-23 [2] CRAN (R 4.3.2)
# globals            0.16.2     2022-11-21 [2] CRAN (R 4.3.2)
# glue               1.7.0      2024-01-09 [2] CRAN (R 4.3.2)
# goftest            1.2-3      2021-10-07 [2] CRAN (R 4.3.2)
# gridExtra          2.3        2017-09-09 [2] CRAN (R 4.3.2)
# gtable             0.3.4      2023-08-21 [2] CRAN (R 4.3.2)
# harmony          * 1.2.0      2023-11-29 [2] CRAN (R 4.3.2)
# here             * 1.0.1      2020-12-13 [2] CRAN (R 4.3.2)
# hms                1.1.3      2023-03-21 [2] CRAN (R 4.3.2)
# htmltools          0.5.7      2023-11-03 [2] CRAN (R 4.3.2)
# htmlwidgets        1.6.4      2023-12-06 [2] CRAN (R 4.3.2)
# httpuv             1.6.14     2024-01-26 [2] CRAN (R 4.3.2)
# httr               1.4.7      2023-08-15 [2] CRAN (R 4.3.2)
# ica                1.0-3      2022-07-08 [2] CRAN (R 4.3.2)
# igraph             2.0.1.9008 2024-02-09 [2] Github (igraph/rigraph@39158c6)
# irlba              2.3.5.1    2022-10-03 [2] CRAN (R 4.3.2)
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
# magrittr           2.0.3      2022-03-30 [2] CRAN (R 4.3.2)
# MASS               7.3-60.0.1 2024-01-13 [3] CRAN (R 4.3.2)
# Matrix             1.6-5      2024-01-11 [3] CRAN (R 4.3.2)
# matrixStats        1.2.0      2023-12-11 [2] CRAN (R 4.3.2)
# mime               0.12       2021-09-28 [2] CRAN (R 4.3.2)
# miniUI             0.1.1.1    2018-05-18 [2] CRAN (R 4.3.2)
# munsell            0.5.0      2018-06-12 [2] CRAN (R 4.3.2)
# nlme               3.1-164    2023-11-27 [3] CRAN (R 4.3.2)
# parallelly         1.36.0     2023-05-26 [2] CRAN (R 4.3.2)
# patchwork        * 1.2.0      2024-01-08 [2] CRAN (R 4.3.2)
# pbapply            1.7-2      2023-06-27 [2] CRAN (R 4.3.2)
# pillar             1.9.0      2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig          2.0.3      2019-09-22 [2] CRAN (R 4.3.2)
# plotly             4.10.4     2024-01-13 [2] CRAN (R 4.3.2)
# plyr               1.8.9      2023-10-02 [2] CRAN (R 4.3.2)
# png                0.1-8      2022-11-29 [2] CRAN (R 4.3.2)
# polyclip           1.10-6     2023-09-27 [2] CRAN (R 4.3.2)
# progressr          0.14.0     2023-08-10 [2] CRAN (R 4.3.2)
# promises           1.2.1      2023-08-10 [2] CRAN (R 4.3.2)
# purrr            * 1.0.2      2023-08-10 [2] CRAN (R 4.3.2)
# R6                 2.5.1      2021-08-19 [2] CRAN (R 4.3.2)
# ragg               1.2.7      2023-12-11 [2] CRAN (R 4.3.2)
# RANN               2.6.1      2019-01-08 [2] CRAN (R 4.3.2)
# RColorBrewer       1.1-3      2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp             * 1.0.12     2024-01-09 [2] CRAN (R 4.3.2)
# RcppAnnoy          0.0.22     2024-01-23 [2] CRAN (R 4.3.2)
# RcppHNSW           0.6.0      2024-02-04 [2] CRAN (R 4.3.2)
# readr            * 2.1.5      2024-01-10 [2] CRAN (R 4.3.2)
# reshape2           1.4.4      2020-04-09 [2] CRAN (R 4.3.2)
# reticulate         1.35.0     2024-01-31 [2] CRAN (R 4.3.2)
# RhpcBLASctl        0.23-42    2023-02-11 [2] CRAN (R 4.3.2)
# rlang              1.1.3      2024-01-10 [2] CRAN (R 4.3.2)
# ROCR               1.0-11     2020-05-02 [2] CRAN (R 4.3.2)
# rprojroot          2.0.4      2023-11-05 [2] CRAN (R 4.3.2)
# RSpectra           0.16-1     2022-04-24 [2] CRAN (R 4.3.2)
# Rtsne              0.17       2023-12-07 [2] CRAN (R 4.3.2)
# scales             1.3.0      2023-11-28 [2] CRAN (R 4.3.2)
# scattermore        1.2        2023-06-12 [2] CRAN (R 4.3.2)
# sctransform        0.4.1      2023-10-19 [2] CRAN (R 4.3.2)
# sessioninfo      * 1.2.2      2021-12-06 [2] CRAN (R 4.3.2)
# Seurat           * 5.0.1      2023-11-17 [2] CRAN (R 4.3.2)
# SeuratObject     * 5.0.1      2023-11-17 [2] CRAN (R 4.3.2)
# shiny              1.8.0      2023-11-17 [2] CRAN (R 4.3.2)
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
# systemfonts        1.0.5      2023-10-09 [2] CRAN (R 4.3.2)
# tensor             1.5        2012-05-05 [2] CRAN (R 4.3.2)
# textshaping        0.3.7      2023-10-09 [2] CRAN (R 4.3.2)
# tibble           * 3.2.1      2023-03-20 [2] CRAN (R 4.3.2)
# tidyr            * 1.3.1      2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect         1.2.0      2022-10-10 [2] CRAN (R 4.3.2)
# tidyverse        * 2.0.0      2023-02-22 [2] CRAN (R 4.3.2)
# timechange         0.3.0      2024-01-18 [2] CRAN (R 4.3.2)
# tzdb               0.4.0      2023-05-12 [2] CRAN (R 4.3.2)
# utf8               1.2.4      2023-10-22 [2] CRAN (R 4.3.2)
# uwot               0.1.16     2023-06-29 [2] CRAN (R 4.3.2)
# vctrs              0.6.5      2023-12-01 [2] CRAN (R 4.3.2)
# viridisLite        0.4.2      2023-05-02 [2] CRAN (R 4.3.2)
# withr              3.0.0      2024-01-16 [2] CRAN (R 4.3.2)
# xtable             1.8-4      2019-04-21 [2] CRAN (R 4.3.2)
# zoo                1.8-12     2023-04-13 [2] CRAN (R 4.3.2)
# 
# [1] /users/csoto/R/4.3.x
# [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library
