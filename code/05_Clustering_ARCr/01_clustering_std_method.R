#########################################################################
##
## Cluster combining multiple parameters (4 algorithms, 2 resolution thr and 3 knn points)
##
## Input:  CellRangerARC_reanalyze QCed norm-data (Seurat) including dim.red and harmonized assays
## Output: Seurat with WNN clusters
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
args = commandArgs(trailingOnly = TRUE)
# clust_method = c(2,3,4), clust_res = c(0.8, 1, 2), clust_knn = c(20, 30, 40)
clust_method <- args[2]
clust_res <- args[4]
clust_knn <- args[6]

# ## For testing use:
# clust_method = 4 Leiden
# clust_res = 2
# clust_knn = 30

message(
  "Processing ",
  clust_method,
  " at res=",
  clust_res,
  " with knn=",
  clust_knn
)

if (length(clust_method) && length(clust_res) && length(clust_knn)) {
  message(
    "\n ====== Processing clustering with method ",
    clust_method,
    " at resolution=",
    clust_res,
    " with k.nn = ",
    clust_knn,
    " ======\n"
  )
} else {
  message("Input arguments missed")
  stop()
}

## Preparing directories
inputDir <- here(
  "processed-data",
  "03_pseudobulking",
  "cellrangerARC_reanalyze"
)
outputRDS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "01_clustering_std_method"
)
outputCVS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "01_clustering_std_method",
  "cvs_files"
)
plotDir <- here("plots", "05_Clustering_ARCr", "01_clustering_std_method")


## Check if processed_data directory exists, if not create it
if (!dir.exists(inputDir)) {
  message("RDS input data missed!")
  stop()
}
if (!dir.exists(outputRDS_Dir)) {
  dir.create(outputRDS_Dir)
}
if (!dir.exists(outputCVS_Dir)) {
  dir.create(outputCVS_Dir)
}
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}


########################    Initials. (1) load data  ########################

# count_mtx_type <- 'norm_counts'
# Seurat_reduction <- 'Harmony'
# minCells <- 1
# if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts'}
# Seurat_base_name <-  paste0(Seurat_base_name, "_Harmony_ARCr_QCed")
# rds_name <- here(inputDir, paste0(Seurat_base_name, ".rds"))

## Read seurat QC'ed with batch correction on both rna and atac assays
SeuratOBJ <- readRDS(here(
  inputDir,
  "seurat.norm_counts_ARCr_harmony_atac_rna_QCed.rds"
))
SeuratOBJ
class(SeuratOBJ[["ATAC"]])
# [1] "ChromatinAssay"
# attr(,"package")
# [1] "Signac"
message('Seurat object loaded!')

# table(SeuratOBJ$orig.ident)
# S03_Hb_r S04_Hb_r S05_Hb_r S06_Hb_r S07_Hb_r S08_Hb_r S09_Hb_r S10_Hb_r
# 4375     6364     1990     7059     7207     4735     5959     5914
# S11_Hb_r S12_Hb_r
# 8054     4045

message("Clustering ", length(Cells(SeuratOBJ)), " cells")


####### Clustering for RNA

## There are two approaches for running the RNA analysis, in this script we use  the (1) Standard seurat workflow

## (1) Standard Seurat workflow

set.seed(03122024)

# Note data set were previously normalized and scaled

DefaultAssay(SeuratOBJ) <- "RNA"

SeuratOBJ.1 <- FindNeighbors(
  SeuratOBJ,
  dims = 1:30,
  reduction = "integrated.harmony"
)

## get name to save the clustering results
clust_name <- case_when(
  clust_method == 1 ~ "C.louvain",
  clust_method == 2 ~ "C.louvainM",
  clust_method == 3 ~ "C.SLM",
  clust_method == 4 ~ "C.leiden"
)


SeuratOBJ.1 <- FindClusters(
  SeuratOBJ.1,
  resolution = as.integer(clust_res),
  method = "igraph", # is deprecated as of Seurat 5.2.0
  random.seed = 03122024,
  algorithm = as.integer(clust_method),
  cluster.name = clust_name
)
# Returns a Seurat where the idents have been updated with new cluster info
# latest clustering results will be stored in object metadata under 'seurat_clusters'.
# Note that 'seurat_clusters' will be overwritten every time FindClusters is run
head(SeuratOBJ.1@meta.data[[clust_name]])
# tail(SeuratOBJ.1[["seurat_clusters"]], n=3)
# tail(SeuratOBJ.1[[clust_name]], n=3)
# pmatch(SeuratOBJ.1[["seurat_clusters"]], SeuratOBJ.1[[clust_name]])

message("\nSeurat Default Clustering:")
df1 <- as.data.frame(table(Idents(SeuratOBJ)))
message("\nSeurat Re-clustering:")
df2 <- as.data.frame(table(Idents(SeuratOBJ.1)))
df_clusters_compare <- merge(df1, df2, by = "Var1", all = TRUE)
df_compare <- df_clusters_compare |> arrange((Var1)) # |> mutate(diff_size = (`Freq.x`-`Freq.y`))

rm("SeuratOBJ")

message("Number of communities: ", length(table(Idents(SeuratOBJ.1))))

print(df_compare, row.names = FALSE)
#      Var1 Freq.x Freq.y diff_size
# 1     0   5539   5300       239
# 2     1   4255   4805      -550
# 3     2   4236   3291       945

## customize umap assay name accordingly with clustering method
umap_rna_name <- paste0("umap.", gsub("\\C.", "", clust_name))
# umap_rna_name <- case_when(
#   clust_method == "1" ~ "umap.lovain",
#   clust_method == "2" ~ "umap.lovainM",
#   clust_method == "3" ~ "umap.SLM",
#   clust_method == "4" ~ "umap.leiden")

message(
  "Clustering method: ",
  clust_name,
  "\n",
  "Clustering resolution: ",
  as.integer(clust_res),
  "\n",
  "Clustering knn: ",
  clust_knn
)

message(
  "RNA SNN ",
  clust_name,
  " method = found ",
  length(unique(SeuratOBJ.1$seurat_clusters)),
  " clusters"
)


########### Next run ATAC analysis. ###########

## ATAC modality was prepared on ~/Hb_multiome/code/03_pseudobulking/08_harmony_CR_ARCr.R

########### Run WNN on integrated.harmony for rna and lsi (1) and lsi-harmony (2) for atac ###########

## ATAC data are normalized and batch corrected
DefaultAssay(SeuratOBJ) <- "ATAC"
clust_name_atac <- paste0(clust_name, "_atac")
# Reductions(SeuratOBJ.1) # "umap.lsi.integrated"
SeuratOBJ.1 <- FindNeighbors(
  SeuratOBJ.1,
  dims = 2:20,
  reduction = "integrated.lsi.harmony"
)
## Build nearest neighbor graph using ATAC (LSI-Harmony)
SeuratOBJ.1 <- FindClusters(
  SeuratOBJ.1,
  resolution = as.integer(clust_res),
  method = "igraph", # is deprecated as of Seurat 5.2.0
  random.seed = 03122024,
  algorithm = as.integer(clust_method),
  cluster.name = clust_name_atac
)

message(
  "ATAC SNN ",
  clust_name,
  " method = found ",
  length(unique(SeuratOBJ.1$seurat_clusters)),
  " clusters"
)

## verify output
colnames(SeuratOBJ.1@meta.data)

## strategy: WNN on integrated.harmony for rna and integrated.lsi.harmony for atac

## Calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering
## We used the reduction with batch corrected data for both rna and atac

DefaultAssay(SeuratOBJ) <- "RNA"

SeuratOBJ.1 <- FindMultiModalNeighbors(
  SeuratOBJ.1,
  k.nn = as.integer(clust_knn),
  reduction.list = list("integrated.harmony", "integrated.lsi.harmony"),
  dims.list = list(1:30, 2:20)
)
# reduction.list = list("pca", "lsi"),
# reduction.list = list("integrated.harmony", "lsi"),

message("WNN done!")
# head(SeuratOBJ.1@neighbors)
# head(SeuratOBJ.1@graphs$wknn)

SeuratOBJ.1 <- RunUMAP(
  SeuratOBJ.1,
  n.neighbors = as.integer(clust_knn), # Default n.neighbors=30
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

SeuratOBJ.1 <- FindClusters(
  SeuratOBJ.1,
  graph.name = "wsnn",
  method = "igraph",
  random.seed = 03122024,
  resolution = as.integer(clust_res),
  algorithm = as.integer(clust_method)
)

# str(SeuratOBJ.1$wnn.umap)
# table(Idents(SeuratOBJ.1))
colnames(SeuratOBJ.1@meta.data)
unique(SeuratOBJ.1$seurat_clusters)

## reduction = "umap.integrated" | var: umap_rna_name
plt1 <- DimPlot(
  SeuratOBJ.1,
  reduction = "umap.integrated",
  group.by = clust_name,
  label = TRUE,
  label.size = 2.5,
  repel = TRUE
) +
  ggtitle(
    "RNA",
    subtitle = paste(
      clust_name,
      " at res=",
      clust_res,
      " (knn=",
      clust_knn,
      "); Clusters=",
      subtitle = length(table(SeuratOBJ.1[[clust_name]]))
    )
  ) &
  NoLegend()
plt2 <- DimPlot(
  SeuratOBJ.1,
  reduction = "umap.lsi.integrated",
  group.by = clust_name_atac,
  label = TRUE,
  label.size = 2.5,
  repel = TRUE
) +
  ggtitle(
    "ATAC",
    subtitle = paste(
      clust_name,
      " at res=",
      clust_res,
      " (knn=",
      clust_knn,
      "); Clusters=",
      subtitle = length(table(SeuratOBJ.1[[clust_name_atac]]))
    )
  ) &
  NoLegend()
plt3 <- DimPlot(
  SeuratOBJ.1,
  reduction = "wnn.umap",
  group.by = "seurat_clusters",
  label = TRUE,
  label.size = 2.5
) +
  ggtitle(
    "WNN",
    subtitle = paste(
      clust_name,
      " at res=",
      clust_res,
      " (knn=",
      clust_knn,
      "); Clusters=",
      subtitle = length(table(SeuratOBJ.1[["seurat_clusters"]]))
    )
  )

pltALL <- plt1 + plt2 + plt3 & theme(plot.title = element_text(hjust = 0.5)) # & NoLegend()
sufix_name <- paste0("k", clust_knn, "_", clust_name, "_lsi_r", clust_res)
ggsave(
  pltALL,
  filename = here(
    plotDir,
    paste0(
      "seurat.norm_counts_CRr_UMAP_WNN_rnaHarm_atacHarm_",
      sufix_name,
      "_v2.png"
    )
  ),
  height = 7,
  width = 20
)


## Save RDS Object and Find DEG in the WNN

clust_knn
if (
  (clust_knn == 30 || clust_knn == 40) &
    (clust_method == 2 || clust_method == 4)
) {
  rds_name <- here(
    outputRDS_Dir,
    paste0("seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_", sufix_name, ".rds")
  )
  saveRDS(SeuratOBJ.1, file = rds_name)

  message(
    "\nSeurat with WNN with rna-harmony and atac-harmony saved: ",
    basename(rds_name)
  )

  ## Find DEG in the integrated Seurat for ALL clusters

  table(SeuratOBJ.1[["seurat_clusters"]])
  all.markers <- FindAllMarkers(object = SeuratOBJ.1)
  head(all.markers, n = 3)

  cvs_file <- paste0(
    "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_",
    sufix_name,
    "_markers.csv"
  )
  cvs_file <- here(outputCVS_Dir, cvs_file)
  write.csv(all.markers, cvs_file)

  message(
    "\nMarkers from WNN with rna-harmony and atac-harmony saved: ",
    basename(cvs_file)
  )
}

# ########### Second strategy: WNN on integrated.harmony for rna and on lsi for atac
#
# SeuratOBJ.2 <- FindMultiModalNeighbors(SeuratOBJ.1,
#                                        k.nn = as.integer(clust_knn),
#                                        reduction.list = list("integrated.harmony", "lsi"),
#                                        dims.list = list(1:30, 2:20))
# message("WNN done!")
# head(SeuratOBJ.2@neighbors)
# rm("SeuratOBJ.1")
# # head(SeuratOBJ.2@graphs$wknn)
#
# SeuratOBJ.2 <- RunUMAP(SeuratOBJ.2,
#                        n.neighbors = as.integer(clust_knn), # Default n.neighbors=30
#                        nn.name = "weighted.nn",
#                        reduction.name = "wnn.umap",
#                        reduction.key = "wnnUMAP_")
#
# SeuratOBJ.2 <- FindClusters(SeuratOBJ.2,
#                             graph.name = "wsnn",
#                             method = "igraph",
#                             random.seed = 03122024,
#                             resolution = as.integer(clust_res),
#                             algorithm = as.integer(clust_method))
#
# #str(SeuratOBJ.2$wnn.umap)
# table(Idents(SeuratOBJ.2))
#
# ## reduction = "umap.integrated" | var: umap_rna_name
# plt1 <- DimPlot(SeuratOBJ.2, reduction = "umap.integrated", group.by = clust_name,
#                 label = TRUE, label.size = 2.5, repel = TRUE) +
#   ggtitle(paste0("RNA (", clust_name, " at res=", clust_res,")")) & NoLegend()
# plt2 <- DimPlot(SeuratOBJ.2, reduction = "umap.lsi.unintegrated", group.by = clust_name,
#                 label = TRUE, label.size = 2.5, repel = TRUE) +
#   ggtitle(paste0("ATAC (LSI)")) & NoLegend()
# plt3 <- DimPlot(SeuratOBJ.2, reduction = "wnn.umap", group.by = clust_name,
#                 label = TRUE, label.size = 2.5) +
#   ggtitle(paste0("WNN (knn=", clust_knn, ", res=", as.character(clust_res), ")"))
#
# pltALL <- plt1 + plt2 + plt3 & theme(plot.title = element_text(hjust = 0.5))
# sufix_name <- paste0("k", clust_knn, "_", clust_name,"_lsi_r", clust_res)
# ggsave(pltALL, filename = here(plotDir, paste0("seurat.norm_counts_CRr_UMAP_WNN_rnaHarm_atacLSI_", sufix_name, ".png")),
#        height = 7, width = 20)
#
#
# ## Find DEG and save Seurat with ONLY ATAC OUTLIER cells (barcodes) to identify cell types later
#
# ## Calculate and save DEG found in the WNN
# if ((clust_knn==30 || clust_knn==40) & (clust_method==2 || clust_method==4)) {
#   rds_name <- here(outputRDS_Dir, paste0("seurat.norm_counts_CRr_WNN_rnaHarm_atacLSI_", sufix_name, ".rds"))
#   saveRDS(SeuratOBJ.2, file = rds_name)
#   message("\nSeurat with WNN with rna-harmony and atac-lsi saved: ", basename(rds_name))
#
#   ## Find DEG in the integrated Seurat for ALL clusters
#
#   table(SeuratOBJ.2[["seurat_clusters"]])
#   all.markers <- FindAllMarkers(object = SeuratOBJ.2)
#   head(all.markers, n=3)
#
#   cvs_file <- paste0("seurat.norm_counts_CRr_WNN_rnaHarm_atacLSI_", sufix_name,"_markers.csv")
#   cvs_file <- here(outputCVS_Dir, cvs_file)
#   write.csv(all.markers, cvs_file)
#
#   message("\nMarkers from WNN with rna-harmony and atac-lsi saved: ", basename(cvs_file))
# }

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
