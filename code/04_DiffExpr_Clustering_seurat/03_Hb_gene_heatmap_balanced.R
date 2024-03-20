
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
library('dplyr')
library('data.table')
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

count_mtx_type <- 'data_counts'
#count_mtx_type <- 'norm_counts' 
Seurat_reduction <- 'CCA'
#Seurat_reduction <- 'Harmony'

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


## Set pre-selected clusters and marker genes to plot

if (Seurat_reduction=='CCA') {
  
  clust_selected <- c(8, 15, 18)
  ## genes selected from Hb_pilot, defined as Hb general, LHb and MHb
  Hb <- c("MMRN1", "GPR151", "POU4F1") # Hb
  LHb <- c("LINC01876") # LHb
  MHb <- c("CD24", "AC004594.1") #MHb
  markers.to.plot <- c(Hb, LHb, MHb)
  
} else {
  
  clust_selected <- c(6, 12, 14, 15)
  ## genes selected from Hb_pilot, defined as Hb general, LHb and MHb
  ## marker genes pre-selected (from Hb_pilot)
  Hb <- c("MMRN1", "GPR151", "POU4F1") # Hb
  LHb <- c("EPHA5", "TLL1") # LHb
  MHb <- c("AC109466.1", "AC008415.1", "GPR149", "GNG8", "NEUROD1", "RASGRP1", "SLC5A7", "CHRNB3", "SCUBE1", "LINC02143", "CD24") #MHb
  markers.to.plot <- c(Hb, LHb, MHb)
  
}  


## Load Seurat integrated with CCA / Harmony
SeuratOBJ <- get_seurat(rds_name)
#SeuratOBJ@reductions
#colnames(SeuratOBJ@meta.data)
head(SeuratOBJ)


## Set specific clusters
SeuratOBJ_Hb_all <- subset(SeuratOBJ, 
                       subset = seurat_clusters %in% clust_selected)

## calculate number of cells by cluster 
clust_tab <- table(SeuratOBJ_Hb_all$seurat_clusters)
# 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
# 0   0   0   0   0   0 830   0   0   0   0   0 180   0 102  99   0 

## get the minimun cluster size
min_clust <- min(clust_tab[clust_tab > 0])

## Set specific marker genes 
markers.to.plot
x_all <- FetchData(object = SeuratOBJ_Hb_all, vars = markers.to.plot)

# ## get some stats over expr. levels
# summary(x_all)
# # Sort the mean expr values 
# x_all %>%
#   summarise(across(all_of(markers.to.plot), ~ mean(.x, na.rm = TRUE))) -> expr_mean
# t(apply(expr_mean,1,sort))

## Select top cells and order them by total expression level 

top_cells <- unlist(lapply(clust_selected, function(current_cluster) {
  x_current <- FetchData(
    object = subset(SeuratOBJ, subset = seurat_clusters %in% current_cluster),
    vars = markers.to.plot
  )
  total_expr <- rowSums(x_current)
  rownames(x_current)[order(total_expr, decreasing = TRUE)][seq_len(min_clust)]
  
}))

p1 <- DoHeatmap(object = SeuratOBJ_Hb_all, 
          cells = top_cells, 
          features=markers.to.plot, label = TRUE, angle=45)

## select cells based in mean expr value for each gene
# LINC02143 <- WhichCells(object = SeuratOBJ_Hb_all, idents = clust_selected, expression = LINC02143 >= 0.1027878)
# SCUBE1  <- WhichCells(object = SeuratOBJ_Hb_all, idents = clust_selected, expression = SCUBE1 >= 0.1103314)
# CHRNB3  <- WhichCells(object = SeuratOBJ_Hb_all, idents = clust_selected, expression = CHRNB3 >= 0.1492207)

# Combine cells
#subset_cells <- c(LINC02143,    SCUBE1,    CHRNB3)

pdf_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_heatmap_balanced.pdf')
pdf_name <- here('plots/04_DiffExpr_Clustering_seurat', pdf_file)
pdf(file = pdf_name)




## Pseudo bulked version

SeuratOBJ_Hb_all_pseudobulked <- AggregateExpression(SeuratOBJ_Hb_all, return.seurat = TRUE, group.by = c("seurat_clusters", "orig.ident"))
table(SeuratOBJ_Hb_all_pseudobulked$orig.ident)
table(SeuratOBJ_Hb_all_pseudobulked$seurat_clusters)

p2 <- DoHeatmap(object = SeuratOBJ_Hb_all_pseudobulked, 
                features=markers.to.plot, label = TRUE, angle=45, group.by = "seurat_clusters",
                size=4) 
p2 + scale_fill_gradientn(limits = c(-2, 2), colours = PurpleAndYellow(), na.value = "white")

pdf_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_heatmap_balanced_pseudobulked.pdf')
pdf_name <- here('plots/04_DiffExpr_Clustering_seurat', pdf_file)
pdf(file = pdf_name)
dev.off()




############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()