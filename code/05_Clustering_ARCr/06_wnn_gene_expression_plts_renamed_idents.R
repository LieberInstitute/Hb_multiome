########################################################################
## Plot GEX on selected WNN clustering results
##
## Authors. CSC 
## Date. Jan 24, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
library("patchwork")
library("purrr")
library("stringr")
library("here")

## input directories

here()

# Check/create directories
inputCVS_Dir <- here("code", "05_Clustering_ARCr")
# inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method") # clusters not labeled
# inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "05_rename_idents") # clusters renamed to identify Hb-clusters
inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "08_wnn_gene_expression_plts_renamed_idents") # clusters renamed for Spatial-Registration on Visium project
plotDir <- here("plots", "05_Clustering_ARCr", "08_wnn_gene_expression_plts_renamed_idents")

## Check directories
if (!dir.exists(plotDir)) {dir.create(plotDir)}


## Load input with RDS wnn to compare

# WNN clustering results of interest
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1"
# Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2.rds"
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"

wnn_file_names_lst <- here(inputRDS_Dir, Seurat_base_name)

# f_plt_violin <- function(seurat_name){
  seurat_name <- wnn_file_names_lst
  SeuratOBJ <- readRDS(seurat_name)
  DefaultAssay(SeuratOBJ) <- "RNA"
  Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+")) 
  # features <- c("POU4F1", "GPR151", "TAC3")
  features <- c("POU4F1", "GPR151")
  
  plt1 <- VlnPlot(object = SeuratOBJ, layer = "data",
                  features = features,
                  pt.size = 0) +
    labs(x = paste0("WNN: ", Seurat_base_name)) &
    theme(text = element_text(size = 8), 
          axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
          plot.title=element_text(hjust=0.5)) 
    
  plt1 <- plt1 + plot_annotation(paste0("WNN: ", Seurat_base_name), 
                                 caption = 'Cell Ranger ARC reanalize',theme=theme(plot.title=element_text(hjust=0.5)))
  tmp_name <- paste0(Seurat_base_name, "_POU4F1_GPR151_VPlot.pdf")
  ggsave(plt1, filename = here(plotDir, tmp_name), height = 4, width = 17)
  
  # return(plt1)
  
# }


## Load the Seurat with WNN clusters

message("Reading WNN to evalute gene expression of `POU4F1` and `GPR151` on: ", str_extract(Seurat_base_name, regex("C\\.\\w+")))

# Vplot_lst <- f_plt_violin(wnn_file_names_lst[1])

# # we have 4 clustering results selected 
# Vplot_lst <- map(wnn_file_names_lst, ~ f_plt_violin(.x))
# length(Vplot_lst)
# # Vplot_lst[[1]]
# 
# pdf_name <-  here(plotDir, "Vplot_WNN.POU4F1_GPR151.pdf")
# pdf(file = pdf_name)
# par(mfrow=c(2,2))
# Vplot_lst[[1]] / Vplot_lst[[2]] / Vplot_lst[[3]] / Vplot_lst[[4]]
# dev.off()

message('\nPlots saved `', plotDir, '`')

## UMAP: Label clusters on a ggplot2-based scatter plot

#plot <- DimPlot(object = sob, 
plt1 <- DimPlot(SeuratOBJ,  label = TRUE, reduction = "wnn.umap", label.size = 3) + NoLegend() +
  labs(title = paste0("Clusters from WNN: ", Seurat_base_name))
tmp_name <- paste0(Seurat_base_name, "_POU4F1_GPR151_DimPlot.pdf")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 6, width = 6)

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


## Single cell heatmap of feature expression
# DoHeatmap(subset(SeuratOBJ, downsample = 100), features = features, size = 3) + labs(title = Seurat_base_name)


## Plot coverage plot for specific genes 
features <- "GPR151"
features <- "POU4F1"
features <- "CDH4"
newname_clusters <- levels(SeuratOBJ)
for (idx in newname_clusters) { if (nchar(idx) <= 4) { no_hb_clust <- append(no_hb_clust, idx) } }
no_hb_clust <- sort(c(unlist(no_hb_clust)))
hb_clusters <- sort(newname_clusters[! newname_clusters %in% c(no_hb_clust)])


plt1 <- CoveragePlot(
  object = SeuratOBJ,
  region = features,
  features = features,
  expression.assay = "RNA",
  extend.upstream = 500,
  extend.downstream = 500,
  idents = hb_clusters
)  +
  labs(title = paste0("Clusters from WNN: ", seurat_name)) +
  theme(text = element_text(size = 8),
        axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
        plot.title=element_text(hjust=0.5))
tmp_name <- paste0(seurat_name, "_", features, "_CoveragePlt.pdf")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 6, width = 6)



## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


