########################################################################
## Plot visualizations for MERGED WNN clusters
##
## Authors. CSC
## Date. Jan 24, 2024
## Recommended resources on interactive mode: srun --pty --mem=80GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("ggplot2")
# library("pheatmap")
# library("bluster")
# library("viridisLite")
# library("patchwork")
# library("ggplotify")
# library("gridExtra")
library("tidyverse")
# library("stringr")
library("here")

## input directories
## use RDS with clusters renamed, also used for Spatial-Registration on Visium project

## input dirs
inputRDS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "05_rename_idents"
)
## output dirs
plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "09_wnn_geneExp_merged_plts_annotated.R"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}

Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
title_name <- str_extract(seurat_name, regex("C\\.\\w*\\_r2"))

# Load Seurat
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
all_clusters <- levels(SeuratOBJ)
all_clusters
#[1] "No-Hb_merged" "MHb_merged"   "LHb_merged" 

# assign new identities to the Seurat object
Idents(SeuratOBJ) <- SeuratOBJ$merged_cluster
levels(SeuratOBJ)
#Levels: No-Hb_merged MHb_merged LHb_merged

table(SeuratOBJ$merged_cluster)
# LHb_merged   MHb_merged No-Hb_merged 
# 6883        10944        37875 


## Plot DimPlot merged clusters

plt1 <- DimPlot(SeuratOBJ, 
                reduction = "wnn.umap",
                group.by = "merged_cluster",
                label = FALSE) +
    labs(title = paste0("WNN clusters merged: ", title_name)) +
    theme(
        text = element_text(size = 8),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5)
    )
f_name <- paste0(Seurat_base_name, "_MERGED_clusters_DimPlotV2.png")
ggsave(plt1, filename = here(plotDir, f_name), height = 6, width = 6)


## Plot DotPlot merged clusters

features <- c("POU4F1", "GPR151", "TAC3")
#title_name <- str_extract(seurat_name, regex("C\\.\\w*\\_r2"))

plt1 <- DotPlot(SeuratOBJ, features = features, 
                group.by = "merged_cluster") + #, dot.scale = 3
    labs(title = paste0("WNN clusters merged: ", title_name)) +
    #theme_minimal() +
    theme(
        text = element_text(size = 12),
        axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5)
    )
#Warning: Scaling data with a low number of groups may produce misleading results
f_name <- paste0(Seurat_base_name, "_MERGED_clusters_DotPlot.pdf")
ggsave(plt1, filename = here(plotDir, f_name), height = 4, width = 6)


message("Plots done!")



## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

