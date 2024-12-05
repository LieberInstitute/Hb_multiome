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
library('ggplot2')
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
plt1 <- ElbowPlot(SeuratOBJ, ndims = 30, reduction = "pca") + ggtitle("PCA reduction")
plt2 <- ElbowPlot(SeuratOBJ, ndims = 30, reduction = "integrated.harmony") + ggtitle("Harmony reduction")
ggsave((plt1 + plt2), filename = here(plotDir, paste0(Seurat_base_name, '_elbow_pca_harmony.png')), height = 5, width = 15) 

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
SeuratOBJ <- RunPCA(SeuratOBJ, features = VariableFeatures(object = SeuratOBJ))
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
df_compare <- df_clusters_compare |> arrange((Var1)) |> mutate(diff_size = (`Freq.x`-`Freq.y`))
print(df_compare, row.names = FALSE)
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
## save as dataframe (raw.names = FALSE)

SeuratOBJ_2 <- RunUMAP(SeuratOBJ_2, dims = 1:10, reduction = "pca", reduction.name = "umap.lovain")

colnames(SeuratOBJ_2@meta.data)
tail(SeuratOBJ_2[["seurat_clusters"]], n=3)
tail(SeuratOBJ_2[["C.Lovain"]], n=3)

plt1 <- DimPlot(SeuratOBJ_2, reduction = "umap.lovain") + labs(title = paste0("UMAP Lovain res=0.8")) +
  labs(subtitle = "CellRangerARC-reanalyze Human Hb")
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

DefaultAssay(SeuratOBJ_2) <- "ATAC"

# Run term frequency inverse document frequency (TF-IDF) normalization on a matrix.
# SeuratOBJ_2 <- RunTFIDF(SeuratOBJ_2,
#                       method = 1,  # computes log(𝑇𝐹×𝐼𝐷𝐹).
#                       scale.factor = 10000)

SeuratOBJ_2 <- FindTopFeatures(SeuratOBJ_2, 
                             min.cutoff = 'q5',
                             verbose = TRUE)
# Set 'q5':  include 95% most common features as the VariableFeatures.
# Set 10:  include features with >10 total counts in the set of VariableFeatures

SeuratOBJ_2 <- RunSVD(SeuratOBJ_2)

Reductions(SeuratOBJ_2)
SeuratOBJ_2 <- RunUMAP(SeuratOBJ_2, 
                     reduction = 'lsi', 
                     dims = 2:50, 
                     reduction.name = "umap.atac", reduction.key = "atacUMAP_")



## Calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering

SeuratOBJ_2 <- FindMultiModalNeighbors(SeuratOBJ_2, 
                                     reduction.list = list("pca", "lsi"),
                                     dims.list = list(1:50, 2:50))

SeuratOBJ_2 <- RunUMAP(SeuratOBJ_2, 
                     nn.name = "weighted.nn", 
                     reduction.name = "wnn.umap", 
                     reduction.key = "wnnUMAP_")

SeuratOBJ_2 <- FindClusters(SeuratOBJ_2, 
                          graph.name = "wsnn", 
                          algorithm = 1, # 1 = original Louvain algorithm
                          verbose = FALSE)

Reductions(SeuratOBJ_2)
table(Idents(SeuratOBJ_2))
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
# 3745 3679 2997 2601 2266 2266 2200 2108 2080 1907 1844 1696 1620 1582 1497 1412
# 16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31
# 1334 1296 1281 1244 1115 1091 1086  983  918  874  874  822  707  690  657  644
# 32   33   34   35   36   37   38   39   40   41   42   43   44
# 619  569  542  508  409  406  269  257  251  246  217  211   82

print(as.data.frame(table(Idents(SeuratOBJ_2))), row.names = FALSE)

plt1 <- DimPlot(SeuratOBJ_2, reduction = "umap.lovain", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
plt2 <- DimPlot(SeuratOBJ_2, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC") 
plt3 <- DimPlot(SeuratOBJ_2, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN") 
pltALL <- plt1 + plt2 + plt3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
ggsave(pltALL, filename = here(plotDir, paste0(Seurat_base_name, '_UMAPS_PCA_LSI_WNN.png')), height = 7, width = 20) 


## Save RDS Object
rds_name <- here(outputRDS_Dir, paste0(Seurat_base_name,'_WNN_clusters.rds'))
saveRDS(SeuratOBJ_2, file = rds_name)
message('Seurat with ATAC clusters saved!')  


## Find DEG in the integrated Seurat for ALL clusters
table(SeuratOBJ_2[["seurat_clusters"]])
all.markers <- FindAllMarkers(object = SeuratOBJ_2)
#head(all.markers, n=3)

cvs_file <- paste0(Seurat_base_name, "_markers_WNN.csv") 
cvs_file <- here(cvsDir, cvs_file)
write.csv(all.markers, cvs_file)

message(" FindAllMarkers done!")

message("All tasks done!")


library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
