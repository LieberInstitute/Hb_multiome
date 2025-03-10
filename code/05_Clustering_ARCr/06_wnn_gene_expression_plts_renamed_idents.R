########################################################################
## Plot GEX on selected WNN clustering results
##
## Authors. CSC 
## Date. Jan 24, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
########################################################################

library("Seurat")
library("ggplot2")
library("patchwork")
library("purrr")
library("stringr")
library("here")

## input directories

# Check/create directories
inputCVS_Dir <- here("code", "05_Clustering_ARCr")
# inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method") # clusters not labeled
inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "05_rename_idents") # clusters renamed to identify Hb-clusters
plotDir <- here("plots", "05_Clustering_ARCr", "04_wnn_gene_expression.R")

## Check directories
if (!dir.exists(plotDir)) {dir.create(plotDir)}
if (!dir.exists(outputCVS_Dir)) {dir.create(outputCVS_Dir)}


## Load input with RDS wnn to compare

# WNN clustering results of interest
# Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1"
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2.rds"

# tmp_dir <- here(inputCVS_Dir, "input_wnn_rds_names.txt")
# wnn_file_names = readLines(tmp_dir)[1:4] # only resolution r1, excluded r2

# any("r1" %in% wnn_file_names)

# wnn_file_names_lst <- here(inputRDS_Dir, paste0(wnn_file_names, ".rds"))
wnn_file_names_lst <- here(inputRDS_Dir, Seurat_base_name)

f_plt_violin <- function(seurat_name){
  # seurat_name = wnn_file_names_lst
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
    
  #plt1 <- plt1 + plot_annotation(paste0("WNN: ", Seurat_base_name), 
  #                               caption = 'Cell Ranger ARC reanalize',theme=theme(plot.title=element_text(hjust=0.5)))
  
  tmp_name <- paste0("Vplot_POU4F1_GPR151_", Seurat_base_name,".pdf")
  ggsave(plt1, filename = here(plotDir, tmp_name), height = 4, width = 17)
  
  return(plt1)
  
}


## Load the Seurat with WNN clusters

message("Reading WNN to evalute gene expression of `POU4F1` and `GPR151` on: ", str_extract(Seurat_base_name, regex("C\\.\\w+")))

# Vplot_lst <- f_plt_violin(wnn_file_names_lst[1])

# we have 4 clustering results selected 
Vplot_lst <- map(wnn_file_names_lst, ~ f_plt_violin(.x))
length(Vplot_lst)
# Vplot_lst[[1]]

pdf_name <-  here(plotDir, "Vplot_WNN.POU4F1_GPR151.pdf")
pdf(file = pdf_name)
par(mfrow=c(2,2))
Vplot_lst[[1]] / Vplot_lst[[2]] / Vplot_lst[[3]] / Vplot_lst[[4]]
dev.off()

message('\nPlots saved `', plotDir, '`')

## UMAP: Label clusters on a ggplot2-based scatter plot
sob <- readRDS(wnn_file_names_lst[1])
DefaultAssay(sob) <- "RNA"
plot <- DimPlot(object = sob, 
                reduction = "wnn.umap")
LabelClusters(plot = plot, id = 'ident') + labs(title = Seurat_base_name)

# Feature plot - visualize feature expression in low-dimensional space
# Calculate feature-specific contrast levels based on quantiles of non-zero expression.
# Particularly useful when plotting multiple markers
Reductions(sob)
FeaturePlot(sob, features = features,
            reduction = "wnn.umap",
            min.cutoff = "q10", max.cutoff = "q90") + labs(title = Seurat_base_name)
            # min.cutoff = 1, max.cutoff = 3)

# Visualize co-expression of two features simultaneously
FeaturePlot(sob, features = features, blend = TRUE) # + labs(title = Seurat_base_name)

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(sob, features = features) + RotatedAxis() + labs(title = Seurat_base_name)

# Single cell heatmap of feature expression
DoHeatmap(subset(sob, downsample = 100), features = features, size = 3) + labs(title = Seurat_base_name)


# # DoHeatmap now shows a grouping bar, splitting the heatmap into groups or clusters. This can
# # be changed with the `group.by` parameter
# DoHeatmap(sob, features = VariableFeatures(sob)[1:100], cells = 1:500, size = 4,
#           angle = 90) + NoLegend()


# # Rename cell identity classes
# # Can provide an arbitrary amount of idents to rename
# levels(pbmc_small)
# #> [1] "0" "1" "2"
# pbmc_small <- RenameIdents(pbmc_small, '0' = 'A', '2' = 'C')
# levels(pbmc_small)
# #> [1] "A" "C" "1"


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


