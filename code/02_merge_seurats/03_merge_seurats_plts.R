########################################################################
## Plot merge Seurat objects contained in an specific directory
##
## Authors. CSC
## Date. Feb 22, 2024
## Last.md: xxx
##
## Input: Merged Seurat RDS Object 
## Output:  Plots before integration for reference    
##
## NOTES: 20G free mem recommended for 20K cells
## For slurm env: sbatch 03_merge_seurats_plts.sh
########################################################################

library('Seurat')              
library('dplyr')
options(tidyverse.quiet = TRUE)
library('tidyverse')
library('ggplot2')
library('patchwork')
library('here')

here::here()


# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/02_merge_seurats/"))) {
    dir.create(here("plots/02_merge_seurats/"))
}



########################    Initials ########################  

## Counts for GEX assay available in Seurat: raw and normalized data
# count_mtx_type_label <- 'data_counts'      
# count_mtx_type_label <- 'norm_counts' 

## get args

args = commandArgs(trailingOnly=TRUE)
## read count mtx type (abs_counts and normalized_counts)
count_mtx_type_label <- args[2]

#message('\nPlotting merged Seurat(s) of the `', count_mtx_type_label, '` layer')

## compose rds base name for merge data
if (count_mtx_type_label=='data_counts') { 
  s_sample <- 'seurat.combined.data_counts' 
  rna_layer <- 'counts'
} else if (count_mtx_type_label=='norm_counts') { # normalized counts
  s_sample <- 'seurat.combined.norm_counts' 
  rna_layer <- 'data'
} else {
  stop()
}

message('\nPlotting `', rna_layer, '` layer\n' )



############ Load Seurat to plot  ############

## read directory with Seurat objects
path_directory_in <- here('processed-data/02_merge_seurats')
rds_path <- paste0(path_directory_in, '/', list.files(path_directory_in, pattern=paste0(s_sample,'.rds'))) #, recursive = TRUE
rds_path

if (length(rds_path) < 1) { stop('\nNone Seurat available.') }

print(rds_path)
SeuratOBJ.combined <-readRDS(rds_path)

print(SeuratOBJ.combined)
DefaultAssay(SeuratOBJ.combined) <- "RNA"

message('\nSeurat combined loaded ', split(table(SeuratOBJ.combined$orig.ident), ','))






############ Plots  UMIs, Genes, ^MT and RIBO levels  ############

message('\nBuilding plots ...')


#colnames(SeuratOBJ.combined@meta.data)

## Violin plot main GEX features
## Note:  VlnPlot is plotting the data from the data slot by default
p1 <- VlnPlot(object = SeuratOBJ.combined, layer = rna_layer, 
              features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
              group.by = 'orig.ident', pt.size = 0) & geom_boxplot() &
  theme(legend.position = 'none',
        axis.text.x = element_text(angle=45, size=8),  
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8)) 

#head( GetAssayData(object = SeuratOBJ.combined[["RNA"]], layer = 'data')[])
#SeuratOBJ.combined@assays$RNA@counts
#head(SeuratOBJ.combined[[]])

## Plot Genes and UMIs by density per cell 
df_genes_per_cell <- as.data.frame(SeuratOBJ.combined[[]])
p2 <- df_genes_per_cell %>%
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(size=8),  
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)) +
  ylab("Log10(UMIs)") +
  xlab("Gene-counts") +
  ggtitle("Genes density by cell") + theme(plot.title = element_text(size = 8))

## Correlation btw Genes/UMIs 
p4 <- df_genes_per_cell %>%
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt, group.by = 'orig.ident')) + # MTRatio
  #    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=MTRatio)) + # MTRatio
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + scale_y_log10() +
  theme_classic() +
  theme(axis.text.x = element_text(size=8),  
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8)) +
  facet_wrap(~orig.ident) +
  ylab("log10(nGenes)") +
  xlab("log10(UMIs)") +
  ggtitle('UMIs per Genes by MT levels') + theme(plot.title = element_text(size = 8))

pdf_file <- paste0(s_sample, '_ALLplots_GEX.pdf')
pdf_name <-  here('plots/02_merge_seurats', pdf_file)  
pdf(file = pdf_name)
# Example file name: seurat.combined.data_counts_ALLplots_GEX.pdf

pAll <- p1 / p2 / p4 +
  labs(caption = paste('Note: Gene expression is from ', rna_layer, ' slot.'))

pAll

dev.off()

message('\nPlots saved in plots/02_merge_seurats as `', pdf_file, '`')


## slurm script reproducibility

slurmjobs::job_loop(
  loops = list(type_mtx = c("data_counts", "norm_counts")),
  name = "03_merge_seurats_plts",
  cores = 2,
  create_shell = TRUE
)











# ############ Run PCA in the merged Seurat for reference
# 
# # NOTE timoast comment: https://github.com/satijalab/seurat/issues/3505
# # I'd suggest doing QC and filtering cells on each object before running the integration.
# # Running NormalizeData on the integrated assay will overwrite the integration results.
# 
# ## split the RNA measurements into two layers one for each sample
# SeuratOBJ.combined[["RNA"]] <- split(SeuratOBJ.combined[["RNA"]], f = SeuratOBJ.combined$orig.ident)
# # Warning: Assay RNA changing from Assay to Assay5
# 
# SeuratOBJ.combined <- NormalizeData(SeuratOBJ.combined)
# 
# ## Identifies features that are outliers on a 'mean variability plot'.
# ## vst method (default): First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess).
# 
# SeuratOBJ.combined <- FindVariableFeatures(SeuratOBJ.combined, selection.method = "vst") 
# # Fits a line to the relationship of log(variance) and log(mean) using local polynomial regression...
# 
# all.genes <- rownames(SeuratOBJ.combined)
# SeuratOBJ.combined <- ScaleData(SeuratOBJ.combined, features = all.genes)
# ## Scale data. Perform “LogNormalize” method to the GEX for each cell by the total expression multiply by a scale factor (10,000 by default), and log-transforms the result
# 
# ## Run a PCA dimensionality reduction
# SeuratOBJ.combined <- RunPCA(SeuratOBJ.combined)
# #SeuratOBJ.combined@reductions
# 
# ## Save more variable features at top 10,20,50 and 100 VF 
# save_VFeatures(SeuratOBJ.combined, 'pca')
# 
# p1 <- ElbowPlot(SeuratOBJ.combined)
# png_file <- paste0(s_sample, '_PCAelbow.png')
# png_name <- here('plots/02_merge_seurats', png_file)  
# ggsave(p1, filename = png_name, height = 4, width = 5)
# 
# 
# 
# 
# ## Cluster the samples before corrections for further reference
# 
# ## Compute nearest neighbor graph + SNN
# 
# SeuratOBJ.combined <- FindNeighbors(SeuratOBJ.combined, dims = 1:30, reduction = "pca")
# # Compute Louvain
# SeuratOBJ.combined <- FindClusters(SeuratOBJ.combined, 
#                                    resolution = 1, cluster.name = "unintegrated_clusters")
# #unique(SeuratOBJ.combined$seurat_clusters)
# #unique(SeuratOBJ.combined$unintegrated_clusters)
# 
# SeuratOBJ.combined <- RunUMAP(SeuratOBJ.combined, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# head(SeuratOBJ.combined, n=2)
# #SeuratOBJ.combined@reductions
# 
# # create and save UMAP-PCA plots grouped by sample and clusters and splitted side-by-side
# plot_clust(SeuratOBJ.combined, paste0(s_sample, '_umap.unintegrated_'), 'umap.unintegrated', 'seurat_clusters')
# 
# # visualize more variable features in a heatmap
# p1 <- DimHeatmap(SeuratOBJ.combined, reduction = 'pca', nfeatures = 30)
# png_file <- paste0(s_sample, '_umap.unintegrated_pca_heatmap.png')
# png_name <- here('plots/02_merge_seurats', png_file)
# ggsave(p1, filename = png_name) #, height = 5, width = 5
# 
# # Save RDS Object
# rds_name <- here('processed-data/02_merge_seurats', paste0(s_sample, '_PCA.rds'))
# # .../seurat.combined.data_counts_PCA.rds
# saveRDS(SeuratOBJ.combined, file = rds_name)
# message('Seurat combined saved in ', rds_name)   
# 





############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()


# > print('Reproducibility information:')
# [1] "Reproducibility information:"
# > # Last modification
#   > Sys.time()
# [1] "2024-04-04 21:42:26 EDT"
# > #"2023-04-04 12:42:26 EDT"
#   > proc.time()
# user   system  elapsed 
# 29.026    2.299 1961.593 
# > options(width = 120)
# > session_info()
# abind              1.4-5      2016-07-21 [2] CRAN (R 4.3.2)
# beeswarm           0.4.0      2021-06-01 [2] CRAN (R 4.3.2)
# cli                3.6.2      2023-12-11 [2] CRAN (R 4.3.2)
# cluster            2.1.6      2023-12-01 [3] CRAN (R 4.3.2)
# codetools          0.2-19     2023-02-01 [3] CRAN (R 4.3.2)
# colorspace         2.1-0      2023-01-23 [2] CRAN (R 4.3.2)
# cowplot            1.1.3      2024-01-22 [2] CRAN (R 4.3.2)
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
# ggbeeswarm         0.7.2      2023-04-29 [2] CRAN (R 4.3.2)
# ggplot2          * 3.4.4      2023-10-12 [2] CRAN (R 4.3.2)
# ggrastr            1.0.2      2023-06-01 [2] CRAN (R 4.3.2)
# ggrepel            0.9.5      2024-01-10 [2] CRAN (R 4.3.2)
# ggridges           0.5.6      2024-01-23 [2] CRAN (R 4.3.2)
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
# mgcv               1.9-1      2023-12-21 [3] CRAN (R 4.3.2)
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
# RANN               2.6.1      2019-01-08 [2] CRAN (R 4.3.2)
# RColorBrewer       1.1-3      2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp               1.0.12     2024-01-09 [2] CRAN (R 4.3.2)
# RcppAnnoy          0.0.22     2024-01-23 [2] CRAN (R 4.3.2)
# RcppHNSW           0.6.0      2024-02-04 [2] CRAN (R 4.3.2)
# readr            * 2.1.5      2024-01-10 [2] CRAN (R 4.3.2)
# reshape2           1.4.4      2020-04-09 [2] CRAN (R 4.3.2)
# reticulate         1.35.0     2024-01-31 [2] CRAN (R 4.3.2)
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
# slurmjobs          1.2.1      2024-03-20 [1] Github (LieberInstitute/slurmjobs@3fbddb2)
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
