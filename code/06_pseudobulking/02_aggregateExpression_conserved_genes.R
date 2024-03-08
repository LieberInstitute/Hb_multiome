########################################################################
## Aggregate expression for RNA assay for CCA and Harmony correction
## Authors. CSC
## Date. March 5th, 2024
## Last.Adaptation: xxx
##
## Input: Seurat integrated object with samples S1 and S2 after CCA correction
## Output:  
##
## NOTES:
## For slurm env: runsrun --x11 --pty --partition=interactive bash
########################################################################

library('Seurat')                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
## Additional required packages for aggregation
library('multtest')
library('metap')
## Additional required packages for customize clusters
library('scCustomize')
library('magrittr')
library('tidyverse')
# for a much faster version to run FindMarkers() install these packages:
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')
library('presto')
## Packages to plot
library('ggplot2')
library('patchwork')
library('cowplot')
theme_set(theme_cowplot())

#library(SeuratDisk)                             
library(here)

here::here()

if (!packageVersion("Seurat")=='4.9.9.9060') {
    stop
    message('This pipeline was implemented with Seurat v5 and Signac v1.11+ ')
    message('You need the laterst Seurat v5 (‘4.9.9.9060’)')
    message('Current available repository on: https://satijalab.org/seurat/articles/install.html  ') }

# Check if processed_data directory exists, if not create it
if (!dir.exists(here("processed-data/06_pseudobulking/"))) {
    dir.create(here("processed-data/06_pseudobulking/"))
}
# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/06_pseudobulking/"))) {
    dir.create(here("plots/06_pseudobulking/"))
}
# Check if directory to store results exists, if not create it
if (!dir.exists(here("processed-data/06_pseudobulking/csv_files"))) {
  dir.create(here("processed-data/06_pseudobulking/csv_files"))
}

#source(here("code/functions_custom", "remote_plot_functions.R"))    # Call to plot GEX assay
#source(here("code/functions_custom", "remote_filtering_functions.R"))   # Call functions to subset the Seurat object

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

# get QC violin plots
plot_violinQC <- function(seuratOBJ, sfeature, stitle) {
    p1 <- VlnPlot(object = seuratOBJ, features = sfeature, 
                  group.by = 'orig.ident', pt.size = 0) & geom_boxplot() &
        theme(legend.position = 'none',
              axis.text.x = element_text(angle=0, hjust=1, size=8),  #10
              axis.text.y = element_text(size=8), 
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) #&
    #labs(title = "", x = 'Samples', y ="")
    ggtitle(stitle)
    return(p1)
}



##### load pre-existing Seurat objects

## Compose Seurat object name
if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.combined.data_counts_PCA' } else { Seurat_base_name <- 'seurat.combined.norm_counts_PCA' }
if (Seurat_reduction=='CCA') {
  rds_name <- here('processed-data/04_merge_seurats', paste0(Seurat_base_name, '_CCA.rds'))
} else {
  rds_name <- here('processed-data/04_merge_seurats', paste0(Seurat_base_name, '_Harmony.rds'))
}
rds_name
# ~/seurat.combined.data_counts_PCA_Harmony.rds
# ~/seurat.combined.data_counts_PCA_CCA.rds

SeuratOBJ <- get_seurat(rds_name)
# var. for testing: SeuratOBJ2 <- SeuratOBJ
# An object of class Seurat 
# 36601 features across 17994 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
# 3 layers present: data, counts, scale.data
# 4 dimensional reductions calculated: pca, umap.unintegrated, integrated.cca, umap

# verification of the integration
table(SeuratOBJ$orig.ident)


# S1_Hb_KDM S2_Hb_KDM 
# 8178      9816 
head(colnames(SeuratOBJ))
tail(colnames(SeuratOBJ))
colnames(SeuratOBJ@meta.data)
# [1] "orig.ident"                 "atac_peak_region_fragments"
# [3] "atac_fragments"             "nCount_RNA"                
# [5] "nFeature_RNA"               "log10GenesPerUMI"          
# [7] "percent.mt"                 "percent.ribo"              
# [9] "MTRatio"                    "unintegrated_clusters"     
# [11] "seurat_clusters"            "RNA_snn_res.1"   



######################  Customize clusters  ######################


oldIdent <- levels(Idents(SeuratOBJ))
# [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17" "18"

newIdent <- paste("C", 0:(length(oldIdent)-1), sep = "_")
# [1] "C_0"  "C_1"  "C_2"  "C_3"  "C_4"  "C_5"  "C_6"  "C_7"  "C_8"  "C_9"  "C_10" "C_11" "C_12" "C_13" "C_14" "C_15"
# [17] "C_16" "C_17" "C_18

# rename clusters to make them more readable
# require scCustomize/Wrapper funtion to rename clusters
SeuratOBJ <- Rename_Clusters(SeuratOBJ, new_idents = newIdent,
                         meta_col_name = "seurat_clusters.renamed")
# count cells by clusters
table(Idents(SeuratOBJ))
# C_0  C_1  C_2  C_3  C_4  C_5  C_6  C_7  C_8  C_9 C_10 C_11 C_12 C_13 C_14 C_15 C_16 C_17 C_18 
# 4122 1625 1615 1330 1311 1250 1100 1026  821  787  683  565  564  515  273  163  136   75   33 
#View(table(Idents(SeuratOBJ)))

head(SeuratOBJ)


# ## Trying to change orig.ident meta.data
# library(stringr)
# 
# ## Confirm how many samples do we have
# # suffixes <- str_extract(string = colnames(SeuratOBJ), pattern = "[:digit:]$")
# # unique(suffixes)
# 
# unique(SeuratOBJ@meta.data$orig.ident)
# 
# # Create dataframe by sample that contains matching orig.ident code
# meta_by_sample <- tibble::tribble(
#   ~orig.ident,  ~sample_name,
#   1, "S1_Hb",
#   2, "S2_Hb" 
# )
# 
# # Change orig.ident column to factor so that it can be joined later
# meta_by_sample$orig.ident <- as.factor(meta_by_sample$orig.ident)
# 
# # Pull existing meta data where samples are specified by orig.ident and remove everything but orig.ident
# OBJ_meta <- SeuratOBJ@meta.data %>% 
#   select(orig.ident) %>% 
#   rownames_to_column("barcodes")
# 
# # Use full join with object meta data in x position so that by sample meta dataframe is propagated across the by cell meta dataframe from the object.  And then remove orig.ident because it's already present in object meta data.
# full_new_meta <- full_join(x = OBJ_meta, y = meta_by_sample) %>% 
#   column_to_rownames("barcodes") %>% 
#   select(-orig.ident)
# 
# # Use AddMetaData to add new meta data to object
# OBJ <- AddMetaData(object = OBJ, metadata = full_new_meta)



######################  Identify conserved cell type markers ######################
## Implementation from: https://satijalab.org/seurat/articles/integration_introduction.html#identify-conserved-cell-type-markers 

## Run in an integrated Seurat
SeuratOBJ[["RNA"]] <- JoinLayers(SeuratOBJ[["RNA"]])

## unique(Idents(SeuratOBJ))
## Hb.markers <- FindConservedMarkers(SeuratOBJ, ident.1 = "Clust_0", grouping.var = "orig.ident", verbose = FALSE)
## head(nk.markers)

# ## Determine the number of clusters
## https://github.com/satijalab/seurat/issues/6076

# num_clusters <- max(as.numeric(as.character(
#   SeuratOBJ@meta.data$seurat_clusters)))
# 
# ## Cycle through each cluster finding the conserved markers**
# ## Store each dataframe of markers in the misc slot**
# for (i in 0:num_clusters) {
#   SeuratOBJ@misc$temp <- FindConservedMarkers(SeuratOBJ, ident.1 = i, grouping.var = "orig.ident", min.cells.group = 0)
#   names(gene.conditions@misc)[names(gene.conditions@misc)=="temp"] <-
#     paste0(names(gene.conditions), ".cluster_", i, ".markers")
# }

## Plot conserved cell type markers with Doplot() 

# unique(Idents(SeuratOBJ))
# markers.to.plot <- c("MMRN1", "HTR2C", "EPHA5", "GPR151", "POU4F1", 
#                      "AC109466.1", "AC008415.1", "GPR149", "GNG8", "LINC01876", "TLL1", "CD24", "AC004594.1")
# DotPlot(SeuratOBJ, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") +
#   RotatedAxis()
# 
# DotPlot(SeuratOBJ, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
#   RotatedAxis()



######################. Identify differential expressed genes across conditions ######################

## We use AggregateExpression() to aggregate cells of a similar type and condition together to create “pseudobulk” profiles

## To avoid issue when having few cells need to adjust the minimum number of cells
## For example, if there is a cluster "15" that has 0 cells, the function will skip that cluster with a warning (that's perfect). Also if the number of cells is between min.cells.groups (default = 3) and 0, an error is thrown and it stops working. That is why I previously remove from the Seurat Object the cells of the clusters with 3 or less cells for each condition/sample. 
few_cells_samples <- unique(SeuratOBJ@meta.data$orig.ident)
few_cells <- vector()

for (i in 1:length(few_cells_samples)) {   # remove cellstype w/ less than 1 cells in each sample/condition
  few_cells_tmp <- table(SeuratOBJ@meta.data$seurat_clusters.renamed[SeuratOBJ@meta.data$orig.ident == few_cells_samples[i]]) <= 1
  few_cells_tmp <- names(few_cells_tmp)[few_cells_tmp == "TRUE"]
  few_cells <- c(few_cells,few_cells_tmp)
}

message(' Clusters with less than 1 cell: ', length(few_cells))
# > few_cells
# [1] "18"

clusters <- sort(unique(SeuratOBJ@meta.data$seurat_clusters.renamed))
`%notin%` <- Negate(`%in%`) 
clusters <- clusters[clusters %notin% few_cells]  # need to check CSC

colnames(SeuratOBJ@meta.data)
SeuratOBJ@assays
# data, counts, scale.data

aggregate_ifnb <- AggregateExpression(SeuratOBJ, 
                                      assays = 'RNA',
                                      #group.by = c("orig.ident", "seurat_clusters.renamed"), 
                                      group.by = c("orig.ident", "seurat_clusters"), 
                                      return.seurat = TRUE)
# Defaults to: normalization.method = "LogNormalize", scale.factor = 10000
# If return.seurat = TRUE, aggregated values are placed in the 'counts' layer of the returned object

#aggregate_ifnb
# An object of class Seurat 
# 36601 features across 37 samples within 1 assay 
# Active assay: RNA (36601 features, 0 variable features)
# 3 layers present: counts, data, scale.data

table(Idents(SeuratOBJ))
# C_0  C_1  C_2  C_3  C_4  C_5  C_6  C_7  C_8  C_9 C_10 C_11 C_12 C_13 C_14 C_15 C_16 C_17 C_18 
# 4122 1625 1615 1330 1311 1250 1100 1026  821  787  683  565  564  515  273  163  136   75   33 

DEG.response <- FindAllMarkers(SeuratOBJ, 
                               test.use = "wilcox", #'t' for Student's t-test
                               verbose = TRUE)
head(DEG.response, n = 5)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster   gene
# NPAS3      0 -1.8361509 0.078 0.736         0     C_0  NPAS3
# QKI        0 -1.1883844 0.176 0.830         0     C_0    QKI
# ZBTB20     0 -1.6780801 0.068 0.710         0     C_0 ZBTB20
# CADM2      0 -0.6540443 0.229 0.847         0     C_0  CADM2
# MAGI2      0 -0.7831094 0.142 0.748         0     C_0  MAGI2


cvs_name <- paste0(Seurat_base_name,'_', Seurat_reduction, '_Allmarkers.csv')
cvs_name <- here('processed-data/06_pseudobulking/csv_files', cvs_name)
write.csv(DEG.response, cvs_name, row.names=FALSE)
# ~/processed-data/05_DiffExpr_Clustering_Seurat/csv_files/seurat.combined.data_counts_PCA_CCA_Allmarkers.csv"

## Save integrated object with DEG calculated
rds_name <- paste0(Seurat_base_name,'_', Seurat_reduction, '_pseudobulk.rds')
rds_name <- here('processed-data/06_pseudobulking', rds_name)
# file name: seurat.combined.data_counts_PCA_CCA_pseudobulk.rds
saveRDS(SeuratOBJ, file = rds_name)
message('Seurat combined saved in ', rds_name)   


# All this chunck moved to next script. CSC

# ######################. Several visualizations  ######################
# 
# markers.to.plot <- c("MMRN1", "HTR2C", "EPHA5", "GPR151", "POU4F1")
# markers.to.plot <- c("AC109466.1", "AC008415.1", "GPR149", "GNG8")
# markers.to.plot <- c("LINC01876", "TLL1", "CD24", "AC004594.1")
# markers.to.plot <- c("HTR2C")
# 
# ## Plot DEG in aggregate data
# 
# unique(Idents(SeuratOBJ))
# DotPlot(SeuratOBJ, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
#   RotatedAxis()
# # DotPlot(SeuratOBJ, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") +
# #   RotatedAxis()
# 
# # FeaturePlot(SeuratOBJ, features = genes.to.label , split.by = "orig.ident", max.cutoff = 3,
# #             cols = c("grey","red"), reduction = "integrated.cca")
# 
# # Run umap
# SeuratOBJ <- RunUMAP(SeuratOBJ, dims = 1:30, reduction = "integrated.cca")
# SeuratOBJ@reductions
# 
# # Plot in umap features for LHb/MHb marker genes
# FeaturePlot(SeuratOB, features = markers.to.plot , split.by = "orig.ident", max.cutoff = 3,
#             cols = c("grey","red"), reduction = "umap")
# 
# # Plot Violin plots for the same LHb/MHb marker genes
# plots <- VlnPlot(SeuratOBJ, features = markers.to.plot, split.by = "orig.ident", group.by = "seurat_clusters",
#                  pt.size = 0, combine = FALSE)
# wrap_plots(plots = plots, ncol = 1)
# 
# 
# DoHeatmap(
#   SeuratOBJ,
#   features = NULL,
#   cells = NULL,
#   group.by = "orig.ident",
#   group.bar = TRUE,
#   group.colors = NULL,
#   disp.min = -2.5,
#   disp.max = NULL,
#   slot = "scale.data",
#   assay = NULL,
#   label = TRUE,
#   size = 5.5,
#   hjust = 0,
#   vjust = 0,
#   angle = 45,
#   raster = TRUE,
#   draw.lines = TRUE,
#   lines.width = NULL,
#   group.bar.height = 0.02,
#   combine = TRUE
# )
# 
# 
# 
# png_file <- paste0(s_sample, '_integrated.cca_pca_heatmap.png')
# png_name <- here('plots/04_merge_seurats', png_file)
# ggsave(p1, filename = png_name, height = 5, width = 10)





# INTEGRATION methods for Seurat V5:  https://satijalab.org/seurat/articles/seurat5_integration (Oct 31, 2023)
# https://satijalab.org/seurat/articles/integration_introduction.html (Nov 16, 2023)



############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()


# # Last modification
# Sys.time()
# #"2023-04-04 12:42:26 EDT"
# proc.time()
# options(width = 120)
# session_info()
# > library("sessioninfo")
# > print('Reproducibility information:')
# [1] "Reproducibility information:"
# > # Last modification
#   > Sys.time()
# [1] "2024-03-08 15:31:50 EST"
# > #"2023-04-04 12:42:26 EDT"
#   > proc.time()
# user   system  elapsed 
# 767.643   98.692 2702.232 
# > options(width = 120)
# > session_info()
# 
# cluster            2.1.6      2023-12-01 [3] CRAN (R 4.3.2)
# codetools          0.2-19     2023-02-01 [3] CRAN (R 4.3.2)
# colorspace         2.1-0      2023-01-23 [2] CRAN (R 4.3.2)
# cowplot          * 1.1.3      2024-01-22 [2] CRAN (R 4.3.2)
# curl               5.2.0      2023-12-08 [2] CRAN (R 4.3.2)
# data.table       * 1.15.0     2024-01-30 [2] CRAN (R 4.3.2)
# deldir             2.0-2      2023-11-23 [2] CRAN (R 4.3.2)
# desc               1.4.3      2023-12-10 [2] CRAN (R 4.3.2)
# devtools           2.4.5      2022-10-11 [1] CRAN (R 4.3.2)
# digest             0.6.34     2024-01-11 [2] CRAN (R 4.3.2)
# dotCall64          1.1-1      2023-11-28 [2] CRAN (R 4.3.2)
# dplyr            * 1.1.4      2023-11-17 [2] CRAN (R 4.3.2)
# ellipsis           0.3.2      2021-04-29 [2] CRAN (R 4.3.2)
# fansi              1.0.6      2023-12-08 [2] CRAN (R 4.3.2)
# fastDummies        1.7.3      2023-07-06 [2] CRAN (R 4.3.2)
# fastmap            1.1.1      2023-02-24 [2] CRAN (R 4.3.2)
# fitdistrplus       1.1-11     2023-04-25 [2] CRAN (R 4.3.2)
# forcats          * 1.0.0      2023-01-29 [2] CRAN (R 4.3.2)
# fs                 1.6.3      2023-07-20 [2] CRAN (R 4.3.2)
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
# later              1.3.2      2023-12-06 [2] CRAN (R 4.3.2)
# lattice            0.22-5     2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval           0.2.2      2019-03-15 [2] CRAN (R 4.3.2)
# leiden             0.4.3.1    2023-11-17 [2] CRAN (R 4.3.2)
# lifecycle          1.0.4      2023-11-07 [2] CRAN (R 4.3.2)
# limma              3.58.1     2023-10-31 [2] Bioconductor
# listenv            0.9.1      2024-01-29 [2] CRAN (R 4.3.2)
# lmtest             0.9-40     2022-03-21 [2] CRAN (R 4.3.2)
# lubridate        * 1.9.3      2023-09-27 [2] CRAN (R 4.3.2)
# magrittr         * 2.0.3      2022-03-30 [2] CRAN (R 4.3.2)
# MASS               7.3-60.0.1 2024-01-13 [3] CRAN (R 4.3.2)
# mathjaxr           1.6-0      2022-02-28 [1] CRAN (R 4.3.2)
# Matrix             1.6-5      2024-01-11 [3] CRAN (R 4.3.2)
# matrixStats        1.2.0      2023-12-11 [2] CRAN (R 4.3.2)
# memoise            2.0.1      2021-11-26 [2] CRAN (R 4.3.2)
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
# patchwork        * 1.2.0      2024-01-08 [2] CRAN (R 4.3.2)
# pbapply            1.7-2      2023-06-27 [2] CRAN (R 4.3.2)
# pillar             1.9.0      2023-03-22 [2] CRAN (R 4.3.2)
# pkgbuild           1.4.3      2023-12-10 [2] CRAN (R 4.3.2)
# pkgconfig          2.0.3      2019-09-22 [2] CRAN (R 4.3.2)
# pkgload            1.3.4      2024-01-16 [2] CRAN (R 4.3.2)
# plotly             4.10.4     2024-01-13 [2] CRAN (R 4.3.2)
# plotrix            3.8-4      2023-11-10 [2] CRAN (R 4.3.2)
# plyr               1.8.9      2023-10-02 [2] CRAN (R 4.3.2)
# png                0.1-8      2022-11-29 [2] CRAN (R 4.3.2)
# polyclip           1.10-6     2023-09-27 [2] CRAN (R 4.3.2)
# presto           * 1.0.0      2024-03-08 [1] Github (immunogenomics/presto@31dc97f)
# processx           3.8.3      2023-12-10 [2] CRAN (R 4.3.2)
# profvis            0.3.8      2023-05-02 [2] CRAN (R 4.3.2)
# progressr          0.14.0     2023-08-10 [2] CRAN (R 4.3.2)
# promises           1.2.1      2023-08-10 [2] CRAN (R 4.3.2)
# ps                 1.7.6      2024-01-18 [2] CRAN (R 4.3.2)
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
# remotes            2.4.2.1    2023-07-18 [2] CRAN (R 4.3.2)
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
# statmod            1.5.0      2023-01-06 [2] CRAN (R 4.3.2)
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
# urlchecker         1.0.1      2021-11-30 [2] CRAN (R 4.3.2)
# usethis            2.2.2      2023-07-06 [2] CRAN (R 4.3.2)
# utf8               1.2.4      2023-10-22 [2] CRAN (R 4.3.2)
# uwot               0.1.16     2023-06-29 [2] CRAN (R 4.3.2)
# vctrs              0.6.5      2023-12-01 [2] CRAN (R 4.3.2)