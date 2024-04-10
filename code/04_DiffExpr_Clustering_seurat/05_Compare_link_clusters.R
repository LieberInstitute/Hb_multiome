
########################################################################
## Compare link clusters and plot in a heatmap
## INPUT: Seurat object with cluster information
## 
## OUPUT: Heatmap with linkage clusters
## 
## Date: Apr 9th, 2024 
########################################################################

# load libraries
library('Seurat')
library("SingleCellExperiment")
library("dynamicTreeCut")
library("here")
#library("ggplot2")
library("bluster")
library("ComplexHeatmap")

library(here)

here::here()



############        Initials      ############


reduction_name <- 'Harmony'

count_mtx_type <- 'data_counts'



## read directory with Seurat objects

if (count_mtx_type=='data_counts') { 
  Seurat_base_name <- 'seurat.combined.data_counts_PCA'
} else { 
  Seurat_base_name <- 'seurat.combined.norm_counts_PCA' 
}

## Compose Seurat object name processed before
#rds_name1 <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_CCA.rds'))
rds_name <- here('processed-data/02_merge_seurats', paste0(Seurat_base_name, '_Harmony.rds'))

## Check if plots directory exists, if not create it

if (!dir.exists(here("plots/04_DiffExpr_Clustering_seurat/"))) {
  dir.create(here("plots/04_DiffExpr_Clustering_seurat/"))
}
rds_path_out <- here('plots/04_DiffExpr_Clustering_seurat')
rds_path_out



############ Load Seurat objects with clustering data  ############


if ( length(rds_name) < 1 ) { stop('\nNone seurats available.') }

message('Losing Seurats with clustering: \n', rds_name1, ' and \n', rds_name2)


#SeuratOBJ_CCA <-readRDS(rds_name1) #CCA
SeuratOBJ <-readRDS(rds_name) # Harmony
#str(SeuratOBJ_CCA)
#head(SeuratOBJ@reductions$pca) 

#if (SeuratOBJ@reductions %in% c('integrated.harmony','b')) { print(SeuratOBJ$integrated.harmony ) }

#head(SeuratOBJ@reductions$pca) 

# active clusters 
unique(Idents(SeuratOBJ))
#Idents(SeuratOBJ_CCA)[1:1000]





############ Link clusters and build a graph-plot if linkage  ############

# Convert Seurat to sce object; as not supported with v5. I must to convert the assays to v3 assays

Seurat_SCE <- SeuratOBJ
Seurat_SCE[["RNA"]] <- as(Seurat_SCE[["RNA"]], Class="Assay")
sce_h <- as.SingleCellExperiment(Seurat_SCE)
#str(sce_h)
#sce_h$seurat_clusters

clusters <- list(
  Seurat0.8 = sce_h$seurat_clusters,
  Seurat1 = sce_h$RNA_snn_res.1
  )

g <- linkClusters(clusters)
plot(g)

igraph::cluster_walktrap(g)

# Results as a matrix, for two clusterings:
linkClustersMatrix(clusters[[1]], clusters[[2]], denominator="union")





############ Link clusters and build a heatmap of clustering linkage  ############


## Compute the jaccard matrices for comparing clusters like at
## https://github.com/LieberInstitute/spatialDLPFC/blob/14a1f253a92e43c01fec3cc3077a9b2cf9ce9fc0/code/analysis/08_spatial_registration/03_jaccard.R#L95-L109


# Show everything in a data frame:
colData(sce_h)[,"seurat_clusters",drop=FALSE]

levels(sce_h$seurat_clusters)
levels(sce_h$RNA_snn_res.1)


## Create a matrix comparing clusters obtained with k=20 vs k=15
sce_h[['seurat_clusters']]
jacc.mat.08vs1 <-
  with(
    colData(sce_h),
    linkClustersMatrix( sort(sce_h[['seurat_clusters']]), sort(sce_h[['RNA_snn_res.1']]) )
  )

#colnames(jacc.mat.08vs1)

#Heatmap(jacc.mat.08vs1)

## Mark 0s as NAs
#jacc.mat.08vs1[jacc.mat.08vs1 == 0] <- NA

# #Colors for resolution 0.8
domain_colors_res0.8 <- Polychrome::palette36.colors(17)

## Rows annotation (res0.8)
row_ha <- rowAnnotation(
  df = data.frame(res0.8 = c(1:17)),
  col = list(res0.8 = setNames(domain_colors_res0.8, 1:17)),
  show_legend = c(FALSE)
)

# Colors for resolution 1
domain_colors_res1 <- Polychrome::palette36.colors(17)

## Columns annotation (res1)
col_ha <- HeatmapAnnotation(
  df = data.frame(res1 = c(1:17)),
  col = list(K15 = setNames(domain_colors_res1, 1:17)),
  annotation_name_side = "left",
  show_legend = c(FALSE)
)


## Create and save a heatmap for the comparation between clusters k 20 vs 15
#png(file = file.path(dir_plt, "Clusters_res0.8_vs_res1.png"),  width =600)
Heatmap(
  jacc.mat.08vs1,
  name = "Correspondence",
  right_annotation = row_ha,
  bottom_annotation = col_ha,
  col = viridisLite::plasma(101),
  #na_col = "black",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = paste(reduction_name, "\nCompare clusters res0.8 vs res1"),
) 
#dev.off()







############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()