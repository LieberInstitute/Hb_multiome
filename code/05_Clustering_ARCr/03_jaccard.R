########################################################################
## Compares GEX clustering results against pre-selected WNN clustering results
## INPUT:
##      (1) First Seurat with WNN to compare
##      (2) Second Seurat with WNN to compare
##
## OUPUT:
##      1) Jaccard Index Heatmap
##
## Authors. CSC 
## Date. Dec 11, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
########################################################################

library("Seurat")
library("Signac")
library("scran")
library("bluster")
library("SingleCellExperiment")
#library("spatialLIBD")
library("dplyr")
library("bluster")
library("pheatmap")
#library("ComplexHeatmap")
library("ggplot2")
library("stringr")
library("here")
library("viridisLite")

## input directories

# Check/create directories
rna_inputRDS_Dir  <- here("processed-data", "03_pseudobulking", "cellranger_count")
inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
outputCVS_Dir <- here("processed-data", "05_Clustering_ARCr")
plotDir <- here("plots", "05_Clustering_ARCr", "03_jaccard")

## Check directories
if (!dir.exists(plotDir)) {dir.create(plotDir)}
if (!dir.exists(outputCVS_Dir)) {dir.create(outputCVS_Dir)}


## Load ONLY RNA data to make comparison

#seurat_RDSname_1 <- rna_inputRDS_Dir
if ( !length(list.files(rna_inputRDS_Dir, pattern = "seurat.norm_counts_Harmony_All.rds")==1) ) { message("Seurat object missed!");  stop() }
seurat_RDSname_1 <- here(rna_inputRDS_Dir, "seurat.norm_counts_Harmony_All.rds")
Seurat_base_name_1 <- "CR_count"

## Load input with RDS wnn to compare

## read input arguments ( name of RDS Seurat file with wnn clustering to parse )
Seurat_base_name_2 <- commandArgs(trailingOnly = TRUE)
## Some WNN clustering results of interest
## For testing:
Seurat_base_name_2 <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.leiden_lsi_r1"
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.leiden_lsi_r2
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r2
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvainM_lsi_r1
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvainM_lsi_r2
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.SLM_lsi_r1
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.SLM_lsi_r2


## Load the Seurats with WNN clusters

message("Reading Seurat(s) to compute Jaccard Index on WNN clusters")

## Main arguments to specify which WNN clustering results to compare

Seurat_base_name = "seurat.norm_counts_Harmony_ARCr_QCed"
# knn = 20
knn = 30
# knn = 40
# methodWNN = "C.leiden_lsi_r"
# methodWNN = "C.SLM_lsi_r"
# methodWNN = "C.louvain_lsi_r"
# methodWNN = "C.louvainM_lsi_r"
# res = 0.8
resolution = 1
# res = 1.5
# res = 2

# Seurat_base_name_1 = paste0(Seurat_base_name, "_WNN_k", knn, "_", methodWNN, resolution)
# seurat_RDSname_1 = paste0(Seurat_base_name, ".rds")


## Load FIRST Seurat with RNA clustering 

SeuratOBJ_1 <- readRDS(seurat_RDSname_1) # 84177 cells
message("GEX clustering loaded!\nCells: ", length(Cells(x = SeuratOBJ_1)))
n_clust <- nrow(unique(SeuratOBJ_1[["seurat_clusters"]]))
message("\nSNN `", Seurat_base_name_1, "` containing ", n_clust, " clusters")
table(SeuratOBJ_1[["seurat_clusters"]])
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 5519 3930 3647 3112 2653 2277 2272 2265 2130 2100 2049 2044 1739 1697 1640 1614 
# 16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31 
# 1518 1500 1497 1392 1268 1165  990  833  833  702  616  559  466  292  284  238 
# 32   33   34   35   36   37 
# 236  202  194  104   77   48


## Prepare arguments for SECOND WNN data clustering 

# Seurat_base_name = "seurat.norm_counts_Harmony_ARCr_QCed"
methodWNN = "C.leiden_lsi_r"

# Seurat_base_name_2 = paste0(Seurat_base_name_2, "_WNN_k", knn, "_", methodWNN, resolution)
seurat_RDSname_2 = paste0(Seurat_base_name_2, ".rds")
if ( !length(list.files(inputRDS_Dir, pattern = seurat_RDSname_2)==1) ) { message("Seurat 2 object missed!");  stop() }

## Load SECOND Seurat with WNN clustering 

seurat_RDSname_2 <- here(inputRDS_Dir, seurat_RDSname_2)
SeuratOBJ_2 <- readRDS(seurat_RDSname_2)
message("WNN clustering loaded!\nCells: ", length(Cells(x = SeuratOBJ_2)))
n_clust <- nrow(unique(SeuratOBJ_2[["seurat_clusters"]]))
message("\nWNN `", Seurat_base_name_2, "` containing ", n_clust, " clusters")
table(SeuratOBJ_2[["seurat_clusters"]])
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 5516 5240 3041 3023 2378 2330 2295 2261 2123 2077 2041 1971 1809 1758 1743 1700 
# 17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
# 1694 1641 1485 1396 1256 1193  968  892  831  706  615  556  450  238  202  195 
# 33 
# 78 

## Some fast checking

Reductions(SeuratOBJ_1) 
Reductions(SeuratOBJ_2) 
colnames(SeuratOBJ_1@meta.data)
colnames(SeuratOBJ_2@meta.data)

message("\nStarting approximate-silhouette for evaluating cluster separation ...")

## Prepare function to
## (1) Plot approximate silhouette for evaluating cluster separation
## (2) Identified and save closest neighboring cluster for each cell in each cluster 

#plot_approxSilhouette <- function(sce, name_reduction, name_method, re, k){
plot_approxSilhouette <- function(SObj, fn){

  # sil.approx <- approxSilhouette(reducedDim(sce, name_reduction), clusters=colData(sce)$seurat_clusters)
  sil.approx <- approxSilhouette(Embeddings(SObj, reduction = "integrated.harmony"), clusters = SObj$seurat_clusters)
  #sil.approx
  # DataFrame with 84177 rows and 3 columns
  # cluster    other      width
  # <factor> <factor>  <numeric>
  # 10C_AAACAGCCAATCATGT-1       5        13  0.2307765
  # 10C_AAACAGCCACTTCACT-1       16       2   0.4228183
  # 10C_AAACAGCCAGGACCTT-1       0        6   0.1671839
  sil.data <- as.data.frame(sil.approx)
  #sil.data$closest <- factor(ifelse(sil.data$width > 0, colData(sce)$seurat_clusters, sil.data$other))
  sil.data$closest <- factor(ifelse(sil.data$width > 0, SObj$seurat_clusters, sil.data$other))
  #sil.data$cluster <- colData(sce)$seurat_clusters
  sil.data$cluster <- SObj$seurat_clusters
  
  ## identified the closest neighboring cluster for each cell in each cluster
  # tbl_aprox_sil <- table(Cluster=colData(sce)$seurat_clusters, sil.data$closest)
  tbl_aprox_sil <- table(Cluster = SObj$seurat_clusters, sil.data$closest)
  cvs_file <- paste0("Silhouette_", fn, ".cvs")
  cvs_file <- here(outputCVS_Dir, cvs_file)
  write.csv(tbl_aprox_sil, cvs_file)
  message("approximate silhouette cvs saved!")
  
  plt1 <- ggplot(sil.data, aes(x=cluster, y=width, colour=closest)) +
    ggbeeswarm::geom_quasirandom(method="smiley") + labs(title = fn)
    #+ labs(subtitle = "CellRangerARC-reanalyze Human Hb")
  # ggsave(plt1, filename = here(plotDir, paste0("Silhouette_", file_name,".png")), height = 6, width = 10)
  
  return(plt1)
  
}

# sce.1 <- as.SingleCellExperiment(SeuratOBJ_1, assay = "RNA")
# reducedDimNames(sce.1)
# # [1] "PCA"                "UMAP.UNINTEGRATED"  "INTEGRATED.HARMONY"
# # [4] "UMAP"   
# spe1.red_name <-  "INTEGRATED.HARMONY"
# rm("SeuratOBJ_1")
# 
# sce.2 <- as.SingleCellExperiment(SeuratOBJ_2, assay = "RNA")
# reducedDimNames(sce.2)
# 
# "WNN" %in% reducedDimNames(sce.2)
# spe2.red_name <-  "UMAP.LEIDEN"
# rm("SeuratOBJ_2")

##  compute an approximate silhouette width for each observation and plot Silhouette plot
# plt_rna <- plot_approxSilhouette(sce.1, spe1.red_name, substring(spe1.red_name, 6, nchar(spe1.red_name)), resolution, knn)
# Embeddings(SeuratOBJ_1, reduction = "integrated.harmony")
# SeuratOBJ_1[["seurat_clusters"]]


file_name <- "RNA_Harmony"
plt_rna <- plot_approxSilhouette(SeuratOBJ_1, "RNA_Harmony")

file_name <- str_extract(Seurat_base_name_2, "ARC+.+")
plt_wnn <- plot_approxSilhouette(SeuratOBJ_2, file_name)

plt1 <- plt_rna / plt_wnn
tmp_name <- paste0("RNA_", file_name)
tmp_name <- paste0("Silhouette_", tmp_name, ".png")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 12, width = 10)


## Comparing the clustering data sets 

message("Starting Jaccard Index Processing ...")

#https://bioconductor.org/books/3.14/OSCA.advanced/clustering-redux.html
# sce.1 <- as.SingleCellExperiment(SeuratOBJ_1, assay = "RNA")
# clust.louvain <- sce.1[["seurat_clusters"]]
# sce.2 <- as.SingleCellExperiment(SeuratOBJ_2, assay = "RNA")
# clust.leiden <- sce.2[["seurat_clusters"]]
# tab <- table(Louvain=clust.louvain, Leiden=clust.leiden)

clust.rna <- SeuratOBJ_1$seurat_clusters
levels(clust.rna)
clust.wnn <- SeuratOBJ_2$seurat_clusters
levels(clust.wnn)
table(RNA=clust.rna, WNN=clust.wnn)
tab <- table(RNA=clust.rna, WNN=clust.wnn)

# Error in base::table(...) : all arguments must have the same length

rownames(tab) <- paste("WNN", rownames(tab))
colnames(tab) <- paste("RNA", colnames(tab))

pheatmap(log10(tab+10), color=viridis::viridis(100), cluster_cols=FALSE, cluster_rows=FALSE)



jacc.mat <- linkClustersMatrix(clust.rna, clust.wnn)
rownames(jacc.mat) <- paste("WNN", rownames(jacc.mat))
colnames(jacc.mat) <- paste("RNA", colnames(jacc.mat))
plt1 <- pheatmap(jacc.mat, color=viridis::viridis(100), cluster_cols=FALSE, cluster_rows=FALSE)

tmp_name <- paste0(substring(spe1.red_name, 6, nchar(spe1.red_name)), "_", substring(spe2.red_name, 6, nchar(spe2.red_name)))
tmp_name <- paste0("Jaccard_", Seurat_base_name_2, ".png")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 10, width = 10)


## Compute the jaccard matrices, just like at
## https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/4b94e5bf1986df546bdb8624769e2ab746c23e70/code/05_explore_sce/06_explore_azimuth_annotations.R#L109



## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


