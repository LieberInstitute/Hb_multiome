########################################################################
## Plot GEX on selected WNN clustering results
##
## Authors. CSC 
## Date. Jan 24, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
library("patchwork")
library("purrr")
library("tidyverse")
library("stringr")
library("here")

## input directories

here()

# Check/create directories
inputCVS_Dir <- here("code", "05_Clustering_ARCr")
## clusters renamed for Spatial-Registration on Visium project
inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "08_wnn_gene_expression_plts_renamed_idents") 
plotDir <- here("plots", "05_Clustering_ARCr", "08_wnn_gene_expression_plts_renamed_idents")

## Check directories
if (!dir.exists(plotDir)) {dir.create(plotDir)}


## Load input with RDS wnn to compare

# WNN clustering results of interest. To plot annotated or not annotatted clusters
# For inputRDS_Dir_not_annotated
# inputRDS_Dir_not_annotated <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1"
# For inputRDS_Dir_annotated
# inputRDS_Dir_annotated <- here("processed-data", "05_Clustering_ARCr", "05_rename_idents")
# Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2.rds"

# For inputRDS_Dir, clusters renamed for Spatial-Registration on Visium project
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"

seurat_name <- here(inputRDS_Dir, Seurat_base_name)

SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
levels(SeuratOBJ)
## Levels should be
# [1] "C.05 DD_LHb" "C.07 DD_MHb" "C.10 DD_MHb" "C.11 DD_MHb" "C.14 DD_MHb"
# [6] "C.16 DD_MHb" "C.18 DD_LHb" "C.23 DD_LHb" "C.24 DD_LHb" "C.30 DD_LHb"
# [11] "C.33 DD_LHb" "C.36 DD_MHb" "C.40 DD_LHb" "C.01"        "C.02"       
# [16] "C.03"        "C.04"        "C.06"        "C.08"        "C.09"       
# [21] "C.12"        "C.13"        "C.15"        "C.17"        "C.19"       
# [26] "C.20"        "C.21"        "C.22"        "C.25"        "C.26"       
# [31] "C.27"        "C.28"        "C.29"        "C.31"        "C.32"       
# [36] "C.34"        "C.35"        "C.37"        "C.38"        "C.39"       
# [41] "C.41"        "C.42" 


DefaultAssay(SeuratOBJ) <- "RNA"
Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+")) 

## Prepare Violin Plot on canonical Hb gene-markers

message("Reading WNN to evalute gene expression of `POU4F1` and `GPR151` on: ", str_extract(Seurat_base_name, regex("C\\.\\w+")))
# features <- c("POU4F1", "GPR151", "TAC3")
features <- c("POU4F1", "GPR151")
  
plt1 <- VlnPlot(
  object = SeuratOBJ,
  layer = "data",
  features = features,
  pt.size = 0
) +
  labs(x = paste0("WNN: ", Seurat_base_name)) &
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(hjust = 0.5)
  )

plt1 <- plt1 +
  plot_annotation(
    paste0("WNN: ", Seurat_base_name),
    caption = 'Cell Ranger ARC reanalize',
    theme = theme(plot.title = element_text(hjust = 0.5))
  )
tmp_name <- paste0(Seurat_base_name, "_POU4F1_GPR151_VPlot.pdf")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 4, width = 17)

message('\nViolin plots saved `', plotDir, '`')

## UMAP: Label clusters on a ggplot2-based scatter plot

# Feature plot - visualize feature expression in low-dimensional space
# Calculate feature-specific contrast levels based on quantiles of non-zero expression.
# Particularly useful when plotting multiple markers
#Reductions(sob)
# FeaturePlot(sob, features = features,
#             reduction = "wnn.umap",
#             min.cutoff = "q10", max.cutoff = "q90") # + labs(title = Seurat_base_name)

## Visualize co-expression of two features simultaneously
plt1 <- FeaturePlot(SeuratOBJ, features = features, reduction = "wnn.umap", blend = TRUE) +
  labs(title = paste0("Clusters from WNN: ", Seurat_base_name)) &
  theme(text = element_text(size = 8), 
        axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
        plot.title=element_text(hjust=0.5)) 
tmp_name <- paste0(Seurat_base_name, "_POU4F1_GPR151_FeaturePlot.pdf")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 3, width = 10) 


## Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
plt1 <- DotPlot(SeuratOBJ, features = c(features, "TAC3")) + RotatedAxis()  +
  labs(title = paste0("Clusters from WNN: ", Seurat_base_name)) &
  theme(text = element_text(size = 8), 
        axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
        plot.title=element_text(hjust=0.5)) 
tmp_name <- paste0(Seurat_base_name, "_POU4F1_GPR151_DotPlot.pdf")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 6, width = 6)



## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


