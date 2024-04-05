
########################################################################
## Build convergence plot
## INPUT: Seurat object with clusters
## 
## OUPUT: Heatmap with clusters based on jaccard index
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



############ Compute the jaccard matrices  ############


## Compute the jaccard matrices for comparing clusters like at
## https://github.com/LieberInstitute/spatialDLPFC/blob/14a1f253a92e43c01fec3cc3077a9b2cf9ce9fc0/code/analysis/08_spatial_registration/03_jaccard.R#L95-L109

#------------------ Compare clusters obtained with k 20 vs 15 -----------------#

## Create a matrix comparing clusters obtained with k=20 vs k=15
jacc.mat.20vs15 <-
  with(
    colData(sce),
    linkClustersMatrix(prelimClust.20, prelimClust.15)
  )
## Mark 0s as NAs
#jacc.mat.20vs15[jacc.mat.20vs15 == 0] <- NA

# #Colors for k20
domain_colors_k20 <- Polychrome::palette36.colors(13)

## Rows annotation (k20)
row_ha <- rowAnnotation(
  df = data.frame(K20 = c(1:13)),
  col = list(k20 = setNames(domain_colors_k20, 1:13)),
  show_legend = c(FALSE)
)

# Colors for k15
domain_colors_k15 <- Polychrome::palette36.colors(10)

## Columns annotation (k15)
col_ha <- HeatmapAnnotation(
  df = data.frame(K15 = c(1:10)),
  col = list(K15 = setNames(domain_colors_k15, 1:10)),
  annotation_name_side = "left",
  show_legend = c(FALSE)
)

## Create and save a heatmap for the comparation between clusters k 20 vs 15
png(file = file.path(dir_plt, "Clusters_k20_vs_k15.png"),  width =600)
Heatmap(
  jacc.mat.20vs15,
  name = "Correspondence",
  right_annotation = row_ha,
  bottom_annotation = col_ha,
  col = viridisLite::plasma(101),
  #na_col = "black",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "Compare clusters k20 vs k15",
) 
dev.off()







############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()