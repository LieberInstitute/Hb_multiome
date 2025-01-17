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
library("igraph")
library("reticulate")
# ```
# reticulate::py_module_available(module='leidenalg')
# reticulate::import('leidenalg')
# 
# ```
library("leiden")
#library("leidenAlg")
library("tidyverse")
library('ggplot2')
library("here")

here::here()

## read input arguments
args = commandArgs(trailingOnly=TRUE)
# clust_method = c(2,3,4), clust_res = c(0.8, 1, 2), clust_knn = c(20, 30, 40)
clust_method <- args[2]
clust_res <- args[4]
clust_knn <- args[6]

# ## For testing use:
# clust_method = 4
# clust_res = 1
# clust_knn = 20

if (length(clust_method) && length(clust_res) &&  length(clust_knn)) {
  message("\n ====== Processing clustering with method ", clust_method, " at resolution=", clust_res," with k.nn = ", clust_knn, " ======\n")
} else {
  message("Input arguments missed")
  stop() 
}

## Preparing directories
inputDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze")
outputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
outputCVS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method", "cvs_files")
plotDir <- here("plots", "05_Clustering_ARCr", "01_clustering_std_method")


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
SeuratOBJ

message('Seurat object loaded!')

table(SeuratOBJ$orig.ident)
# S03_Hb_r S04_Hb_r S05_Hb_r S06_Hb_r S07_Hb_r S08_Hb_r S09_Hb_r S10_Hb_r
# 4375     6364     1990     7059     7207     4735     5959     5914
# S11_Hb_r S12_Hb_r
# 8054     4045

message("Clustering ", length(Cells(SeuratOBJ)), " cells")

## Plot before re-cluster data for comparison

Reductions(SeuratOBJ)
plt_elbow1 <- ElbowPlot(SeuratOBJ, ndims = 50, reduction = "pca") + ggtitle("PCA reduction") +
  geom_vline(xintercept = 30, color="red", linetype="dashed")
plt_elbow2 <- ElbowPlot(SeuratOBJ, ndims = 50, reduction = "integrated.harmony") + ggtitle("Harmony reduction") +
  geom_vline(xintercept = 30, color="red", linetype="dashed")
ggsave((plt_elbow1 + plt_elbow2), filename = here(plotDir, paste0(Seurat_base_name, '_rna_elbow.png')), height = 6, width = 15)

plt1 <- DimHeatmap(SeuratOBJ, reduction = 'pca', nfeatures = 30, fast = FALSE) + labs(title = paste0("Heatmap PCA")) +
  labs(subtitle = "CellRangerARC-reanalyze Human Habenula dataset")
plt2 <- DimHeatmap(SeuratOBJ, reduction = "integrated.harmony", nfeatures = 30, fast = FALSE) + labs(title = paste0("Heatmap PCA.Harmony")) +
  labs(subtitle = "CellRangerARC-reanalyze Human Habenula dataset") #
ggsave((plt1/plt2), filename = here(plotDir, paste0(Seurat_base_name, '_pca_harmony_heatmap.png')), height = 8, width = 10)

plt1 <- DimPlot(SeuratOBJ, reduction = "pca") + labs(title = paste0("Clustering of ", length(Cells(SeuratOBJ)), " cells")) +
  labs(subtitle = "CellRangerARC-reanalyze Human Habenula dataset") +
  DimPlot(SeuratOBJ, reduction = "umap.unintegrated") +
  DimPlot(SeuratOBJ, reduction = "integrated.harmony")
ggsave(plt1, filename = here(plotDir, paste0(Seurat_base_name, '_redDim_umaps.png')), height = 6, width = 20)



####### Clustering for RNA

## There are two approaches for running the RNA analysis, in this script we use 
##    (1) Standard seurat workflow
##    (2) SCTransform

## SCTransform model normalization is followed by PCA and UMAP dimensionality reduction
##       Function replaces NormalizeData(), ScaleData(), and FindVariableFeatures()


## (1) Standard Seurat workflow

set.seed(03122024)

# Note data set were previously normalized and scaled 

DefaultAssay(SeuratOBJ) <- "RNA"

## Cluster the cells using original Lovain algorithm

SeuratOBJ_2 <- FindNeighbors(SeuratOBJ, dims = 1:30, reduction = "integrated.harmony")

## get proper name to save the clustering results
clust_name <- case_when(
  clust_method == 1 ~ "C.louvain",
  clust_method == 2 ~ "C.louvainM",
  clust_method == 3 ~ "C.SLM",
  clust_method == 4 ~ "C.leiden")

## Seurat methods:
# 1 = original Louvain algorithm
# 2 = Louvain algorithm with multilevel refinement
# 3 = SLM algorithm
# 4 = Leiden algorithm


SeuratOBJ_2 <- FindClusters(SeuratOBJ_2,
                          resolution = as.integer(clust_res),
                          method = "igraph",
                          random.seed = 06122024,
                          algorithm = as.integer(clust_method),
                          cluster.name = clust_name)
# Returns a Seurat where the idents have been updated with new cluster info
# latest clustering results will be stored in object metadata under 'seurat_clusters'.
# Note that 'seurat_clusters' will be overwritten every time FindClusters is run

message("\nClustering after remove atypical cells")
df1 <- as.data.frame(table(Idents(SeuratOBJ)))
message("\nRe-clustering after remove atypical cells")
df2 <- as.data.frame(table(Idents(SeuratOBJ_2)))
df_clusters_compare <- merge(df1, df2, by = "Var1", all = TRUE)
df_compare <- df_clusters_compare |> arrange((Var1)) # |> mutate(diff_size = (`Freq.x`-`Freq.y`))

message("Number of communities: ", length(table(Idents(SeuratOBJ_2))))

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

## get proper name to save umap accordingly with clustering method
umap_rna_name <- case_when(
  clust_method == "1" ~ "umap.lovain",
  clust_method == "2" ~ "umap.lovainM",
  clust_method == "3" ~ "umap.SLM",
  clust_method == "4" ~ "umap.leiden")

message("Clustering method: ", clust_name, "\n", 
        "Clustering resolution: ", as.integer(clust_res), "\n",
        "Clustering knn: ", clust_knn)

#SeuratOBJ_2 <- RunUMAP(SeuratOBJ_2, dims = 1:20, reduction = "pca", reduction.name = "umap.lovain")
SeuratOBJ_2 <- RunUMAP(SeuratOBJ_2, dims = 1:30, 
                       reduction = "integrated.harmony", 
                       reduction.name = umap_rna_name) # umap.lovain

#colnames(SeuratOBJ_2@meta.data)
tail(SeuratOBJ_2[["seurat_clusters"]], n=3)
tail(SeuratOBJ_2[[clust_name]], n=3)

Reductions(SeuratOBJ_2)
plt1 <- DimPlot(SeuratOBJ_2, reduction = umap_rna_name) + labs(title = paste0(umap_rna_name, " at resolution = ", clust_res)) +
  labs(subtitle = "CellRangerARC-reanalyze Human Hb")
ggsave(plt1, filename = here(plotDir, paste0(Seurat_base_name, "_", umap_rna_name, "_r", clust_res,".png")), height = 8, width = 8)



########### Next run ATAC analysis. ###########

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth

DefaultAssay(SeuratOBJ_2) <- "ATAC"

# Run term frequency inverse document frequency (TF-IDF) normalization on a matrix.
SeuratOBJ_2 <- RunTFIDF(SeuratOBJ_2,
                      method = 1,  # computes log(𝑇𝐹×𝐼𝐷𝐹).
                      scale.factor = 10000)

SeuratOBJ_2 <- FindTopFeatures(SeuratOBJ_2,
                             min.cutoff = 'q5',
                             verbose = TRUE)
# Set 'q5':  include 95% most common features as the VariableFeatures.
# Set 10:  include features with >10 total counts in the set of VariableFeatures

SeuratOBJ_2 <- RunSVD(SeuratOBJ_2)

Reductions(SeuratOBJ_2)
plt_elbow3 <- ElbowPlot(SeuratOBJ_2, ndims = 50, reduction = "lsi") + ggtitle("LSI reduction") +
  geom_vline(xintercept = 20, color="red", linetype="dashed")
ggsave((plt_elbow1 + plt_elbow2 + plt_elbow3), filename = here(plotDir, paste0(Seurat_base_name, '_all_elbows.png')), height = 5, width = 20)

SeuratOBJ_2 <- RunUMAP(SeuratOBJ_2,
                     reduction = 'lsi',
                     dims = 2:20,
                     reduction.name = "umap.atac", reduction.key = "atacUMAP_")

plt1 <- DimPlot(SeuratOBJ_2, reduction = "umap.atac") + labs(title = paste0("ATAC-UMAP (LSI)")) +
  labs(subtitle = "CellRangerARC-reanalyze Human Hb")
ggsave(plt1, filename = here(plotDir, paste0(Seurat_base_name, "_umap.lsi.png")), height = 8, width = 8)



## Calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering

SeuratOBJ_2 <- FindMultiModalNeighbors(SeuratOBJ_2,
                                       #reduction.list = list("pca", "lsi"),
                                       k.nn = as.integer(clust_knn),
                                       reduction.list = list("integrated.harmony", "lsi"),
                                       dims.list = list(1:30, 2:20))

SeuratOBJ_2 <- RunUMAP(SeuratOBJ_2,
                       n.neighbors = as.integer(clust_knn), # Default n.neighbors=30
                       nn.name = "weighted.nn",
                       reduction.name = "wnn.umap",
                       reduction.key = "wnnUMAP_")

SeuratOBJ_2 <- FindClusters(SeuratOBJ_2,
                          graph.name = "wsnn",
                          method = "igraph",
                          random.seed = 06122024,
                          resolution = as.integer(clust_res),
                          algorithm = as.integer(clust_method))

Reductions(SeuratOBJ_2)
#str(SeuratOBJ_2$wnn.umap)
table(Idents(SeuratOBJ_2))

# print(as.data.frame(table(Idents(SeuratOBJ_2))), row.names = FALSE)

plt1 <- DimPlot(SeuratOBJ_2, reduction = umap_rna_name, label = TRUE, label.size = 2.5, repel = TRUE) + 
  ggtitle(paste0("RNA (", clust_name, " at res=", clust_res,")")) & NoLegend() # reduction = "umap.lovain"
plt2 <- DimPlot(SeuratOBJ_2, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + 
  ggtitle(paste0("ATAC (LSI)")) & NoLegend()
plt3 <- DimPlot(SeuratOBJ_2, reduction = "wnn.umap", label = TRUE, label.size = 2.5) + 
  ggtitle(paste0("WNN (knn=", clust_knn, ", res=", as.character(clust_res), ")"))
pltALL <- plt1 + plt2 + plt3 & theme(plot.title = element_text(hjust = 0.5)) # & NoLegend()
#pltALL
sufix_name <- paste0("k", clust_knn, "_", clust_name,"_lsi_r", clust_res)
ggsave(pltALL, filename = here(plotDir, paste0(Seurat_base_name, "_UMAP_WNN_", sufix_name, ".png")), height = 7, width = 20)


## Find DEG and save Seurat with ONLY ATAC OUTLIER cells (barcodes) to identify cell types later

## Save RDS Object
rds_name <- here(outputRDS_Dir, paste0(Seurat_base_name, "_WNN_", sufix_name, ".rds"))
saveRDS(SeuratOBJ_2, file = rds_name)

message("\nSeurat with WNN saved: ", basename(rds_name))
message("WNN completed!")

## Find DEG in the integrated Seurat for ALL clusters

DefaultAssay(SeuratOBJ_2) <- "RNA"

table(SeuratOBJ_2[["seurat_clusters"]])
all.markers <- FindAllMarkers(object = SeuratOBJ_2)
head(all.markers, n=3)

cvs_file <- paste0(Seurat_base_name, "_WNN_", sufix_name,"_markers.csv")
cvs_file <- here(outputCVS_Dir, cvs_file)
write.csv(all.markers, cvs_file)

message("\nMarkers saved: ", basename(cvs_file))

message("\nFindAllMarkers completed!")

message("\nAll tasks done!")


# library("slurmjobs")
# slurmjobs::job_loop(
#   loops = list(clust_method = c("1","2","3","4"),
#                clust_res = c("0.8", "1", "1.5","2"), knn = c("20", "30", "40")),
#   name = "01_clustering_std_method",
#   cores = 2,
#   create_shell = TRUE
# )



library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
