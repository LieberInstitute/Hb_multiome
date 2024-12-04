#########################################################################
##
## Clustering in CellRangerARC_reanalyze filtered data ()
##
## Input:  Seurat with ATAC data (rds/h5), meta-data and fragments
## Output: Seurat objects with clustering with standard Seurat workflow --PCA/LSI/WNN 
## 
## Note (1) I apply standard Seurat workflow, data were not normalize with SCT
## Note (2) For +60k cells request 60G free-mem
##
## Authors. CSC 
## Ref: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis#wnn-analysis-of-10x-multiome-rna-atac
## Copied from my own adaptation: https://github.com/cyntsc/RStatClub_Seurat_Signac/blob/main/code/02_Clustering.R  
## Date: Dec 2024
########################################################################

library("Seurat")                                
library("Signac")   
library("tidyverse")
library("here")

here::here()

## Preparing directories
inputDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze")
outputRDS_Dir <- here("processed-data", "05_Clustering_ARCr")
outputCVS_Dir <- here("processed-data", "05_Clustering_ARCr", "cvs_files")
plotDir <- here("plots", "05_Clustering_ARCr")


## Check if processed_data directory exists, if not create it
if (!dir.exists(inputDir)) { message("RDS input data missed!"); stop() }
if (!dir.exists(outputRDS_Dir)) { dir.create(outputRDS_Dir) }
if (!dir.exists(outputCVS_Dir)) { dir.create(outputCVS_Dir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }


########################    Initials. (1) load data  ########################  

count_mtx_type <- 'norm_counts' 
Seurat_reduction <- 'Harmony' 
minCells <- 1
if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts'}
Seurat_base_name <-  paste0(Seurat_base_name, "_Harmony_ARCr_QCed")
rds_name <- here(inputDir, paste0(Seurat_base_name, ".rds"))

## load seurat QC'ed
SeuratOBJ <- readRDS(rds_name)

message('Seurat object loaded!')   

print(SeuratOBJ)
table(SeuratOBJ$orig.ident)
# S03_Hb_r S04_Hb_r S05_Hb_r S06_Hb_r S07_Hb_r S08_Hb_r S09_Hb_r S10_Hb_r 
# 4375     6364     1990     7059     7207     4735     5959     5914 
# S11_Hb_r S12_Hb_r 
# 8054     4045

message("Clustering ", length(Cells(SeuratOBJ)), " cells")
SeuratOBJ

## Plot before re-cluster data for comparison

Reductions(SeuratOBJ)
#ElbowPlot(SeuratOBJ)
plt1 <- DimHeatmap(SeuratOBJ, reduction = 'pca', nfeatures = 30, fast = FALSE) + labs(title = paste0("Heatmap PCA")) +
  labs(subtitle = "CellRangerARC-reanalyze Human Habenula dataset") #
ggsave(plt1, filename = here(plotDir, paste0(Seurat_base_name, '_pca_heatmap.png')), height = 8, width = 8) 

plt1 <- DimPlot(SeuratOBJ, reduction = "pca") + labs(title = paste0("Clustering of ", length(Cells(SeuratOBJ)), " cells")) +
  labs(subtitle = "CellRangerARC-reanalyze Human Habenula dataset") +
  DimPlot(SeuratOBJ, reduction = "umap.unintegrated") + 
  DimPlot(SeuratOBJ, reduction = "integrated.harmony")
ggsave(plt1, filename = here(plotDir, paste0(Seurat_base_name, '_pca_umaps.png')), height = 6, width = 20) 



######## Clustering for RNA

## There are two approaches for running the RNA analysis
##    (1) Standard seurat workflow
##    (2) SCTransform 

## SCTransform model normalization is followed by PCA and UMAP dimensionality reduction
##       Function replaces NormalizeData(), ScaleData(), and FindVariableFeatures()


## (1) Standard Seurat workflow

set.seed(03122024)

DefaultAssay(SeuratOBJ) <- "RNA"

# SeuratOBJ <- NormalizeData(SeuratOBJ, 
#                            normalization.method = "LogNormalize", 
#                            scale.factor = 10000)
# # In the LogNormalize method, Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor
# 
# tail(SeuratOBJ[["RNA"]]$data, n=3)
# 
# SeuratOBJ <- FindVariableFeatures(SeuratOBJ, 
#                                   selection.method = "vst", 
#                                   nfeatures = 2000) 
# # nfeatures define the top variable features to use
# # only used when selection.method is set to 'dispersion' or 'vst'
# 
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(SeuratOBJ), 10)
# 
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(SeuratOBJ, raster=FALSE)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

## Scales and centers features in the dataset
## If variables are provided in vars.to.regress, they are individually regressed against each feature

# all.genes <- rownames(SeuratOBJ)
# SeuratOBJ <- ScaleData(SeuratOBJ, features = all.genes,
#                        vars.to.regress = NULL)

## Perform linear dimensional reduction

# # We perform PCA on the scaled data
# 
# SeuratOBJ <- RunPCA(SeuratOBJ, features = VariableFeatures(object = SeuratOBJ))
# DimHeatmap(SeuratOBJ, dims = 1, cells = 300, balanced = TRUE)
# VizDimLoadings(SeuratOBJ, dims = 1:2, reduction = "pca")
# DimPlot(SeuratOBJ)  + NoLegend()
# ElbowPlot(SeuratOBJ)


## Cluster the cells using original Lovain algorithm 

SeuratOBJ_2 <- FindNeighbors(SeuratOBJ, dims = 1:10)

SeuratOBJ_2 <- FindClusters(SeuratOBJ_2, 
                          resolution = 0.8,
                          algorithm = 1,
                          cluster.name ='C.Lovain')
# Returns a Seurat where the idents have been updated with new cluster info
# latest clustering results will be stored in object metadata under 'seurat_clusters'. 
# Note that 'seurat_clusters' will be overwritten every time FindClusters is run

message("\nClustering after remove atypical cells")
df1 <- as.data.frame(table(Idents(SeuratOBJ)))
message("\nRe-clustering after remove atypical cells")
df2 <- as.data.frame(table(Idents(SeuratOBJ_2)))
df_clusters_compare <- merge(df1, df2, by = "Var1")
df_clusters_compare |>
  arrange((Var1)) |>
  mutate(diff_size = (`Freq.x`-`Freq.y`))
#      Var1 Freq.x Freq.y diff_size
# 1     0   5539   5300       239
# 2     1   4255   4805      -550
# 3     2   4236   3291       945
# 4     3   3800   3180       620
# 5     4   3607   3014       593
# 6     5   3384   2874       510
# 7     6   3255   2674       581
# 8     7   2440   2570      -130
# 9     8   2825   2562       263
# 10    9   2828   2541       287
# ...
# 30   29    101     93         8
# 31   30     61     80       -19

SeuratOBJ_2 <- RunUMAP(SeuratOBJ_2, dims = 1:10, reduction = "pca", reduction.name = "umap.lovain")

colnames(SeuratOBJ_2@meta.data)
tail(SeuratOBJ_2[["seurat_clusters"]], n=3)
tail(SeuratOBJ_2[["C.Lovain"]], n=3)

plt1 <- DimPlot(SeuratOBJ_2, reduction = "umap.lovain") + labs(title = paste0("UMAP Lovain res=0.8")) +
  labs(subtitle = "CellRangerARC-reanalyze Human Hb") #
ggsave(plt1, filename = here(plotDir, paste0(Seurat_base_name, '_UMAP_lovain_0.8.png')), height = 8, width = 8) 

 
# SeuratOBJ <- FindClusters(SeuratOBJ_2, 
#                           resolution = 2,
#                           algorithm = 1,
#                           cluster.name ='C.Lovain.r2')
# 
# SeuratOBJ <- RunUMAP(SeuratOBJ, dims = 1:10, reduction = "pca", reduction.name = "umap.lovain.2")

# ## run the non-linear dimensional reduction (UMAP/tSNE)
#
# DimPlot(SeuratOBJ_2, reduction = "umap.lovain.2")




########### Next run ATAC analysis. ###########

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth

DefaultAssay(SeuratOBJ) <- "ATAC"

SeuratOBJ <- RunTFIDF(SeuratOBJ,
                      method = 1,
                      scale.factor = 10000)

SeuratOBJ <- FindTopFeatures(SeuratOBJ, 
                             min.cutoff = 'q5',
                             verbose = TRUE)
# Set 'q5':  include 95% most common features as the VariableFeatures.
# Set 10:  include features with >10 total counts in the set of VariableFeatures

SeuratOBJ <- RunSVD(SeuratOBJ)

SeuratOBJ <- RunUMAP(SeuratOBJ, 
                     reduction = 'lsi', 
                     dims = 2:50, 
                     reduction.name = "umap.atac", reduction.key = "atacUMAP_")



## Calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering

SeuratOBJ <- FindMultiModalNeighbors(SeuratOBJ, 
                                     reduction.list = list("pca", "lsi"),
                                     dims.list = list(1:50, 2:50))

SeuratOBJ <- RunUMAP(SeuratOBJ, 
                     nn.name = "weighted.nn", 
                     reduction.name = "wnn.umap", 
                     reduction.key = "wnnUMAP_")

SeuratOBJ <- FindClusters(SeuratOBJ, 
                          graph.name = "wsnn", 
                          algorithm = 3, verbose = FALSE)

SeuratOBJ@reductions

library('ggplot2')

p1 <- DimPlot(SeuratOBJ, reduction = "umap.lovain.2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(SeuratOBJ, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(SeuratOBJ, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
pALL <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

table(Idents(SeuratOBJ))

## Base name to save plots
base_name <- levels(SeuratOBJ$`orig.ident`[1])

png_file <- paste0(base_name,'_CLuster_ALL_mt.png')
png_name <- here('plots/GEX_ATAC_preprocessing', png_file)
ggsave(pALL, filename = png_name, height = 4, width = 8)


## perform sub-clustering on a specific clusterto find additional structure

SeuratOBJ <- FindSubCluster(SeuratOBJ, 
    cluster = 0,
    graph.name = 'wsnn', #Name of graph to use for the clustering algorithm
    subcluster.name = "sub.cluster", #name of sub cluster added in the meta.data
    resolution = 0.5,
    algorithm = 1) #1 = original Louvain algorithm

Idents(SeuratOBJ) <- "sub.cluster"

colnames(SeuratOBJ@meta.data)

table(SeuratOBJ@meta.data[['seurat_clusters']])

table(SeuratOBJ@meta.data[['sub.cluster']])


# add annotations
SeuratOBJ <- RenameIdents(SeuratOBJ, '0_0' = 'CD14-A', '0_1' ='CD14-B')
SeuratOBJ$celltype <- Idents(SeuratOBJ)


p1 <- DimPlot(SeuratOBJ, reduction = "umap.lovain.2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(SeuratOBJ, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(SeuratOBJ, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
pALL2 <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

table(Idents(SeuratOBJ))

## Base name to save plots
base_name <- levels(SeuratOBJ$`orig.ident`[1])

png_file <- paste0(base_name,'_CLuster_ALL2_mt.png')
png_name <- here('plots/GEX_ATAC_preprocessing', png_file)
ggsave(pALL, filename = png_name, height = 4, width = 8)


## Save RDS Object
rds_name <- here('processed-data/GEX_ATAC_preprocessing', paste0(base_name,'_Clusters_mt.rds'))
saveRDS(SeuratOBJ, file = rds_name)
message('Seurat with ATAC clusters saved!')  





library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()
