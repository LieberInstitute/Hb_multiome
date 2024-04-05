
########################################################################
## Build convergence plot
## INPUT: Seurat object with clusters
## 
## OUPUT: Heatmap with clusters based in jacquard index
##
## Date: March 29th, 2024 
########################################################################

# load libraries
library('Seurat')                             
library("here")
library("bluster")
library("ComplexHeatmap")

library(here)

here::here()



############        Initials      ############

## read directory with Seurat objects
if (count_mtx_type=='data_counts') { 
  Seurat_base_name <- 'seurat.combined.data_counts_PCA'
} else { 
  Seurat_base_name <- 'seurat.combined.norm_counts_PCA' 
}
## Compose Seurat object name processed before
rds_name1 <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_CCA.rds'))
rds_name2 <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_Harmony.rds'))
rds_name1
rds_name2

# Check if plots directory exists, if not create it
if (!dir.exists(here("plots/04_DiffExpr_Clustering_seurat/"))) {
  dir.create(here("plots/04_DiffExpr_Clustering_seurat/"))
}
rds_path_out <- here('plots/04_DiffExpr_Clustering_seurat')
rds_path_out



############ Load Seurat to plot  ############


if (length(rds_name1) < 1 || length(rds_name2) < 1) { stop('\nNone seurats available.') }

message('Losing Seurats with clustering: \n', rds_name1, ' and \n', rds_name2)

SeuratOBJ1 <-readRDS(rds_name1)
print(SeuratOBJ1@reductions$integrated.cca)
SeuratOBJ2 <-readRDS(rds_name2)
print(SeuratOBJ2@reductions$integrated.harmony)







############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()