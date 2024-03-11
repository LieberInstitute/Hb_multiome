
########################################################################
## Plot gene expression based on the top20 DEG that match with our Habenula custom marker gene list
## INPUT:
##      1) Seurat object with CCA or Harmony reduction with DEG calculated before pseudobulk by orig.ident and cluster
## 
## OUPUT:
##      1) DotPlot
##
## Authors. CSC 
## Date. March 11, 2024 
########################################################################

# load libraries
library('Seurat')                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
library('tidyverse')
#library('dplyr')
#library('data.table')
#library('magrittr')

library(here)

here::here()

# Check if processed_data directory exists, if not create it
if (!dir.exists(here("processed-data/04_DiffExpr_Clustering_seurat/"))) {
  dir.create(here("processed-data/04_DiffExpr_Clustering_seurat/"))
}

# contains the different marker list 
#source(here("code/functions_custom", "remote_DGE_marker_gene_lists.R"))       # Call functions to read paths

get_seurat <- function(name) {
  
  #Ex. /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/02_merge_seurats/seurat.combined.data_counts_PCA.rds
  sobj <- readRDS(name)
  # verification of the integration
  print(table(sobj$orig.ident))
  # S1_Hb_KDM S2_Hb_KDM
  # 8178      9816
  print(head(sobj, n=2))
  return(sobj)
  
}

#############################           Initials        ################################
############################# Pickup a Marker gene list ################################

s_sample <- ''
count_mtx_type <- 'data_counts'
#count_mtx_type <- 'norm_counts' 
#Seurat_reduction <- 'CCA'
Seurat_reduction <- 'Harmony'

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.combined.data_counts_PCA' } else { Seurat_base_name <- 'seurat.combined.norm_counts_PCA' }

## Compose Seurat object name processed before
if (Seurat_reduction=='CCA') {
  rds_name <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_CCA.rds'))
} else {
  rds_name <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_Harmony.rds'))
}
rds_name
# Ex. ~/Hb_multiome/processed-data/02_merge_seurats/seurat.combined.data_counts_PCA_CCA.rds"
# Ex. ~/seurat.combined.data_counts_PCA_Harmony.rds

## Load Seurat Integrated with cluster information
SeuratOBJ <- get_seurat(rds_name)
SeuratOBJ@reductions

## extract cluster data
#md <- SeuratOBJ@meta.data %>% as.data.table

## extract cluster data
colnames(SeuratOBJ@meta.data)
md <- SeuratOBJ@meta.data %>% as.data.table
## Apply vertical format to unique cluster with number of umis, arrenged by sample and cluster number
mdT <- md[, .N, by = c("orig.ident", "seurat_clusters")] %>%
  arrange(., orig.ident, seurat_clusters, .by_group = FALSE)
df_mdT <- as.data.frame(mdT)


# ######################. Several visualizations  ######################

# Adjusted clusters to display according with the clusters where we find Habneula marker genes
#.        DEG found in the top20 

########## For CCA marker genes
## Subset clusters with Habenula information from the Seurat obj
SeuratOBJ_Hb <- subset(SeuratOBJ, subset=seurat_clusters %in% c(8,15,18))
unique(SeuratOBJ_Hb$seurat_clusters)

## Genes arrenged from Hb general, LHb and MHb
markers.to.plot <- c("MMRN1", # Hb
                     "AC004594.1", "LINC01876", # Lb/Mh
                     "CD24") # MHb

########## For Harmony marker genes
SeuratOBJ_Hb <- subset(SeuratOBJ, subset=seurat_clusters %in% c(6, 12, 14, 15))
unique(SeuratOBJ_Hb$seurat_clusters)

## Genes arrenged from Hb general, LHb and MHb
markers.to.plot <- c("MMRN1", "GPR151", "POU4F1", # Hb
                     "EPHA5", "TLL1", # LHb
                     "AC109466.1", "AC008415.1", "GPR149", "GNG8", "NEUROD1", 
                     "RASGRP1", "SLC5A7", "CHRNB3", "SCUBE1", "LINC02143", "CD24") # MHb


## Plot DEG in aggregate data

unique(Idents(SeuratOBJ_Hb))
p1 <- DotPlot(SeuratOBJ_Hb, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()
# DotPlot(SeuratOBJ_Hb, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") +
#   RotatedAxis()

png_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_DoPlot.png')
png_name <- here('plots/04_DiffExpr_Clustering_seurat', png_file)
ggsave(p1, filename = png_name, height = 5, width = 10)


SeuratOBJ_Hb@reductions
genes.to.label <- c("MMRN1", "GPR151", "POU4F1")
FeaturePlot(SeuratOBJ_Hb, features = markers.to.plot, max.cutoff = 3,   # split.by = "orig.ident",
            cols = c("grey","red"), reduction = "integrated.harmony")




# # Run umap
# SeuratOBJ <- RunUMAP(SeuratOBJ, dims = 1:30, reduction = "integrated.cca")
# SeuratOBJ@reductions
# 
# # Plot in umap features for LHb/MHb marker genes
# FeaturePlot(SeuratOB, features = markers.to.plot , split.by = "orig.ident", max.cutoff = 3,
#             cols = c("grey","red"), reduction = "umap")
# 
# Plot Violin plots for the same LHb/MHb marker genes

library('ggplot2')
library('patchwork')
library('cowplot')
theme_set(theme_cowplot())

p1 <- VlnPlot(SeuratOBJ, features = markers.to.plot, split.by = "orig.ident", group.by = "seurat_clusters",
                 pt.size = 0, combine = FALSE)
p1 <- wrap_plots(plots = plots, ncol = 1)

png_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_VPlots.png')
png_name <- here('plots/04_DiffExpr_Clustering_seurat', png_file)
ggsave(p1, filename = png_name, height = 5, width = 10)

DoHeatmap(
  SeuratOBJ,
  features = markers.to.plot,
  cells = NULL,
  group.by = "orig.ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  vjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)


str(SeuratOBJ)

p1 <- DoHeatmap(
  SeuratOBJ,
  features = markers.to.plot,
  cells = NULL,
  group.by = "seurat_clusters",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",  #"scale.data"
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  vjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)

png_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_DoHeatmap.png')
png_name <- here('plots/04_DiffExpr_Clustering_seurat', png_file)
ggsave(p1, filename = png_name, height = 7, width = 12)

############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()