########################################################################
## Rename WNN Clusters Hb clusters (Idents) with cell-types identified with the annotation
## INPUT:
##      (1) GEX DEG annotation
##      (2) Seurat with wnn 
## OUPUT:
##      (1) Seurat with re-named idents
##      (2) Plots for exploration     
## Authors. CSC 
## Date. Jan, 2024
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
########################################################################

library("Seurat")
library("Signac")
library("purrr")
library("ggplot2")
library("dplyr")
library("stringr")
library("here")

## input directories

inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
inputCVS_Dir_Ann <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze_v3", "cvs_files_markers")
outputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "05_rename_idents")
plotDir <- here("plots", "05_Clustering_ARCr", "05_rename_idents")

## Check directories
if (!dir.exists(plotDir)) {dir.create(plotDir)}
if (!dir.exists(outputRDS_Dir)) {dir.create(outputRDS_Dir)}

## read input arguments ( name of RDS Seurat file with wnn clustering to parse )
Seurat_base_name <- commandArgs(trailingOnly = TRUE)

Seurat_base_name <- Seurat_base_name[[2]]

message("Processing ", Seurat_base_name)


##### (1) Load Seurat with WNN idents given by default 

## For testing: 
# Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r2"

seurat_RDSname <- here(inputRDS_Dir, paste0(Seurat_base_name, ".rds"))
     
SeuratOBJ <- readRDS(seurat_RDSname)

total_cells <- length(Cells(x = SeuratOBJ))

message("Renaming ", nrow(unique(SeuratOBJ[["seurat_clusters"]])), " clusters for ", Seurat_base_name)



##### (2) Load DEG with ident annotation

## extract a shorter name to save files 
tmp_wd <- str_extract(Seurat_base_name, pattern = "WNN\\w*\\.\\w*")
# old name: deg_file <- paste0("ARCr_QCed_", tmp_wd, "_DEG_Top50_WNN_DD_LB_matching_markers_integrated.csv")
deg_file <- paste0(tmp_wd, "_cellTypes_integrated_top50.csv")
# WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r1_cellTypes_integrated_top50.csv
file_name <- here(inputCVS_Dir_Ann, deg_file)
file_ann <- read.csv(file_name)
head(file_ann)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster       gene cell_type
# 1     0  -3.230185 0.222 0.693         0       0      CPNE4      <NA>
# 2     0  -3.687100 0.184 0.640         0       0       VAV3      <NA>
# 3     0  -4.194794 0.029 0.460         0       0 AC119673.2      <NA>


## Prepare a list with cell-types to rename idents (takes as input the annotated cell-types)

## get unique annotated cell-types (clusters) containing the 'Hb' word
get_hb_words <- as.data.frame(file_ann) |> filter(grepl('Hb', cell_type))
get_hb_words <- unique(get_hb_words[c("cluster", "cell_type")])


## Calculate % of cells by cluster

f_get_percent_label <- function(cl) {
  ## for testing:  cl = 0
  clust_size <- length(SeuratOBJ$seurat_clusters[SeuratOBJ[["seurat_clusters"]]==cl])
  ## longer label
  # clust_label <- paste0("[", clustSize = clust_size, " / ", paste0(round(((clust_size*100) / total_cells), 2), "%"), "]")
  ## shorter label
  clust_label <- paste0("(", paste0(round(((clust_size*100) / total_cells), 2), "%"), ")")
  return(clust_label)               
}

v_unique_hb_clust <- unlist(unique(get_hb_words["cluster"]))
cells_percents <- map(v_unique_hb_clust, ~ f_get_percent_label(.x))
df_hb <- data.frame(
  cluster = c(v_unique_hb_clust),
  percent = c(unlist(cells_percents)))
df_hb
# cluster percent
# cluster1        0 (9.91%)
# cluster2        1 (7.06%)
# cluster3        2 (6.55%)


## collapse the redundancy clusters

collapsed_cell_types <- aggregate(cell_type ~ cluster, data = get_hb_words, FUN = function(x) paste(x, collapse = ", "))
## set pad of 2 digits to clusters
collapsed_cell_types$cluster <- paste0("C.", sprintf('%02d', collapsed_cell_types$cluster))
## shorter the cell_type label to 8 characters
v_hb_short_cell_types <- substr(collapsed_cell_types$cell_type, 1, 8)
collapsed_cell_types$cell_type <- c(v_hb_short_cell_types)
collapsed_cell_types$cell_type <- trimws(gsub(",", "", collapsed_cell_types$cell_type))

collapsed_cell_types$cell_type <- paste(collapsed_cell_types$cell_type, df_hb$percent)

head(collapsed_cell_types)
# cluster  cell_type
# 1    C.00     DD_MHb
# 2    C.01     DD_LHb
# 3    C.02     DD_MHb
# 4    C.03     DD_LHb
# 5    C.06 LB_Hb neur
# 6    C.08 LB_Hb neur


##### (3) extract original idents (clusters) and prepare the new ident names to rename Seurat clusters

Seurat_clusterIDS <- as.integer(levels(SeuratOBJ$seurat_clusters))
Seurat_clusterIDS <- paste0("C.", sprintf('%02d',Seurat_clusterIDS))
Seurat_clusterIDS <- as.data.frame(Seurat_clusterIDS)
colnames(Seurat_clusterIDS) = "cluster"

## reassign cluster name according with cell-annotation
Seurat_clusterIDS_new <- merge(Seurat_clusterIDS, collapsed_cell_types, all.x = TRUE)
# cluster  cell_type
# 1     C.00     DD_MHb
# 2     C.01     DD_LHb
# 3     C.02     DD_MHb
# 4     C.03     DD_LHb
# 5     C.04       <NA>
# 6     C.05       <NA>
# 7     C.06 LB_Hb neur
# 8     C.07       <NA>

## replace NAs   
Seurat_clusterIDS_new[is.na(Seurat_clusterIDS_new)] <- " "
head(Seurat_clusterIDS_new)
# cluster      cell_type
# 1    C.00 DD_MHb (9.91%)
# 2    C.01 DD_LHb (7.06%)
# 3    C.02 DD_MHb (6.55%)
# 4    C.03 DD_LHb (5.59%)
# 5    C.04               
# 6    C.05       

v_new_clusterIDS <- trimws(paste(Seurat_clusterIDS_new$cluster, Seurat_clusterIDS_new$cell_type))

levels(SeuratOBJ)
# [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14"
# [16] "15" "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29"
# [31] "30" "31" "32" "33" "34" "35" "36" "37"

new_names <- c(v_new_clusterIDS)
names(new_names) <- levels(SeuratOBJ)

## rename idents 
# Idents(SeuratOBJ)
SeuratOBJ <- RenameIdents(object = SeuratOBJ, new_names)
levels(SeuratOBJ)
# [1] "C.00 DD_MHb (9.91%)"   "C.01 DD_LHb (7.06%)"   "C.02 DD_MHb (6.55%)"  
# [4] "C.03 DD_LHb (5.59%)"   "C.04"                  "C.05"                 
# [7] "C.06 DD_MHb (4.08%)"   "C.07"                  "C.08 DD_MHb (3.82%)" 
#SeuratOBJ@meta.data

message("Clusters renaming done!")

## save RDS
rds_file_name <- here(outputRDS_Dir, paste0(Seurat_base_name, ".rds"))
saveRDS(SeuratOBJ, rds_file_name)

message("New seurat with clusters renamed saved!")


##### (4) Some visualizations

message("Building some plots ...")

# Reductions(SeuratOBJ)
## extract suffix name to give unique name to plots
seurat_name <- str_extract(Seurat_base_name, pattern = "k[3:4]0\\_C\\.\\w*")
# seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r2

plt1 <- DimPlot(SeuratOBJ, label = TRUE, reduction = "wnn.umap", label.size = 3) + NoLegend() +
  labs(title = paste0("**Clusters from WNN: ", seurat_name))
tmp_name <- paste0(seurat_name, "_DimPlot_renamed.pdf")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 6, width = 6)

## Violin plots for the features selected

DefaultAssay(SeuratOBJ) <- "RNA"

features <- c("POU4F1", "GPR151") # "TAC3"
  
plt1 <- VlnPlot(object = SeuratOBJ, layer = "data",
                features = features,
                pt.size = 0) +
  labs(x = paste0("**Clusters from WNN: ", seurat_name)) &
  theme(text = element_text(size = 8), 
        axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7)) 

tmp_name <- paste0(seurat_name, "_POU4F1_GPR151_Violin_plot.pdf")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 4, width = 17)
  

## Feature plot - visualize feature expression in low-dimensional space
Reductions(SeuratOBJ)
# FeaturePlot(SeuratOBJ, features = features, reduction = "wnn.umap")
# Visualize co-expression of two features simultaneously
plt1 <- FeaturePlot(SeuratOBJ, features = features, reduction = "wnn.umap", blend = TRUE) +
  labs(title = paste0("**Clusters from WNN: ", seurat_name)) &
  theme(text = element_text(size = 8), 
        axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
        plot.title=element_text(hjust=0.5)) 
tmp_name <- paste0(seurat_name, "_POU4F1_GPR151_FeaturePlot.pdf")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 3, width = 10)


## Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
plt1 <- DotPlot(SeuratOBJ, features = c(features, "TAC3")) + RotatedAxis()  +
  labs(title = paste0("**Clusters from WNN: ", seurat_name)) &
  theme(text = element_text(size = 8), 
        axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
        plot.title=element_text(hjust=0.5)) 
tmp_name <- paste0(seurat_name, "_POU4F1_GPR151_DotPlot.pdf")
ggsave(plt1, filename = here(plotDir, tmp_name), height = 6, width = 6)

message("Plots completed!")

message("Process completed!")


# library("slurmjobs")
# job_loop(
#   loops = list(clustering_name = c("x1", "x2", "x3", "x4")),
#   name = "05_rename_idents",
#   cores = 2,
#   create_shell = TRUE,
#   partition = "katun"
# )


## Additional plots prepared to  TLDR slides 2025

# features <- "GPR151"
# features <- "POU4F1"
# features <- "CDH4"
# idents_to_plt <- c("C.00 DD_MHb", "C.01 DD_LHb", "C.02 DD_MHb","C.03 DD_LHb","C.04", "C.05", "C.06 DD_MHb", "C.07", "C.08 DD_MHb", "C.09 LB_Hb ne", "C.10 DD_LHb", "C.11 DD_MHb")  
# idents_to_plt <- c("C.00 DD_MHb", "C.01 DD_LHb", "C.02 DD_MHb","C.03 DD_LHb","C.08 DD_MHb","C.28","C.20","C.21","C.22","C.23","C.24","C.25")  
# 
# plt1 <- CoveragePlot(
#   object = SeuratOBJ,
#   region = features,
#   features = features,
#   expression.assay = "RNA",
#   extend.upstream = 500,
#   extend.downstream = 500,
#   idents = idents_to_plt
# )  +
#   labs(title = paste0("**Clusters from WNN: ", seurat_name)) +
#   theme(text = element_text(size = 8), 
#         axis.text.x= element_text(size = 7), axis.text.y= element_text(size = 7),
#         plot.title=element_text(hjust=0.5)) 
# tmp_name <- paste0(seurat_name, "_", features, "_CoveragePlt.pdf")
# ggsave(plt1, filename = here(plotDir, tmp_name), height = 6, width = 6)

