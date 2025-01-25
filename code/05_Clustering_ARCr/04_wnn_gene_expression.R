########################################################################
## Plot GEX on selected WNN clustering results
##
## Authors. CSC 
## Date. Jan 24, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
########################################################################

library("Seurat")
library("ggplot2")
library("stringr")
library("here")

## input directories

# Check/create directories
inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
outputCVS_Dir <- here("processed-data", "05_Clustering_ARCr", "04_wnn_gene_expression.R")
plotDir <- here("plots", "05_Clustering_ARCr", "04_wnn_gene_expression.R")

## Check directories
if (!dir.exists(plotDir)) {dir.create(plotDir)}
if (!dir.exists(outputCVS_Dir)) {dir.create(outputCVS_Dir)}


## Load input with RDS wnn to compare

## read input arguments ( name of RDS Seurat file with wnn clustering to parse )
Seurat_base_name <- commandArgs(trailingOnly = TRUE)
## Some WNN clustering results of interest
## For testing:
Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1"
# seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvainM_lsi_r1"

## Prepare RDS seurat names and file names 
seurat_RDSname <- here(inputRDS_Dir, paste0(Seurat_base_name, ".rds"))
Seurat_base_name <- str_extract(Seurat_base_name, regex("C\\.\\w+")) #C.louvain_lsi_r1

## Load the Seurat with WNN clusters

message("Reading Seurat(s) to evalute gene expression: ",Seurat_base_name)

SeuratOBJ <- readRDS(seurat_RDSname)
Reductions(SeuratOBJ)
message("WNN clustering loaded!\nCells: ", length(Cells(x = SeuratOBJ)))
n_clust <- nrow(unique(SeuratOBJ[["seurat_clusters"]]))
message("\nWNN `", basename(Seurat_base_name), "` containing ", n_clust, " clusters")


## Plot cells expressing POU4F1 and GPR151 on each cluster
library(patchwork)
DefaultAssay(SeuratOBJ) <- "RNA"
features <- c("POU4F1", "GPR151")

plt1 <- VlnPlot(object = SeuratOBJ, layer = "data",
                features = features,
                pt.size = 0) 
plt1 <- plt1 + plot_annotation(paste0("WNN: ", Seurat_base_name), 
                               caption = 'Cell Ranger ARC reanalize',theme=theme(plot.title=element_text(hjust=0.5)))

tmp_name <- paste0("Vplot_WNN.", Seurat_base_name, ".png")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 4, width = 16)



## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


