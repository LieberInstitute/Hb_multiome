library("Seurat")
library("Signac")
library("ggplot2")
library("purrr")
library("dplyr")
library("stringr")
library("here")

## input directories

# Check/create directories
inputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method")
inputCVS_Dir_Ann <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze_v2")
outputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "05_rename_idents")
plotDir <- here("plots", "05_Clustering_ARCr", "05_rename_idents")

## Check directories
if (!dir.exists(plotDir)) {dir.create(plotDir)}
if (!dir.exists(outputRDS_Dir)) {dir.create(outputRDS_Dir)}

## read input arguments ( name of RDS Seurat file with wnn clustering to parse )
# Seurat_base_name <- commandArgs(trailingOnly = TRUE)


##### (1) Load Seurat with WNN idents given by default 

## For testing: Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1"

seurat_RDSname <- here(inputRDS_Dir, paste0(Seurat_base_name, ".rds"))
     
SeuratOBJ <- readRDS(seurat_RDSname)
## verification
length(Cells(x = SeuratOBJ))

message("Renaming ", nrow(unique(SeuratOBJ[["seurat_clusters"]])), " clusters for ", Seurat_base_name)



##### (1) Load DEG with ident annotation

tmp_wd <- str_extract(Seurat_base_name, pattern = "WNN\\w*\\.\\w*")
# WNN_k30_C.louvain_lsi_r1
deg_file <- paste0("ARCr_QCed_", tmp_wd, "_DEG_Top50_WNN_DD_LB_matching_markers_integrated.csv")
# "ARCr_QCed_WNN_k30_C.louvain_lsi_r1_DEG_Top50_WNN_DD_LB_matching_markers_integrated.csv"
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
## collapse the redundancy clusters
collapsed_cell_types <- aggregate(cell_type ~ cluster, data = get_hb_words, FUN = function(x) paste(x, collapse = ", "))
## set pad of 2 digits to clusters
collapsed_cell_types$cluster <- paste0("C.", sprintf('%02d', collapsed_cell_types$cluster))
## shorter the cell_type label to 8 characters
v_hb_short_cell_types <- substr(collapsed_cell_types$cell_type, 1, 8)
collapsed_cell_types$cell_type <- c(v_hb_short_cell_types)
collapsed_cell_types$cell_type <- trimws(gsub(",", "", collapsed_cell_types$cell_type))

head(collapsed_cell_types)
# cluster  cell_type
# 1    C.00     DD_MHb
# 2    C.01     DD_LHb
# 3    C.02     DD_MHb
# 4    C.03     DD_LHb
# 5    C.06 LB_Hb neur
# 6    C.08 LB_Hb neur


##### (3) extract idents (clusters) to prepare the new ident names

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
# cluster cell_type
# 1    C.00    DD_MHb
# 2    C.01    DD_LHb
# 3    C.02    DD_MHb
# 4    C.03    DD_LHb
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
SeuratOBJ <- RenameIdents(object = SeuratOBJ, new_names)
levels(SeuratOBJ)
# [1] "C.00 DD_MHb"   "C.01 DD_LHb"   "C.02 DD_MHb"   "C.03 DD_LHb"  
# [5] "C.04"          "C.05"          "C.06 DD_MHb"   "C.07"         
# [9] "C.08 DD_MHb"   "C.09 LB_Hb ne" "C.10 DD_LHb"   "C.11 DD_MHb"  


## save RDS
rds_file_name <- here(outputRDS_Dir, paste0(Seurat_base_name, ".rds"))
saveRDS(SeuratOBJ, rds_file_name)




##### (4) Some visualizations

# DimPlot(SeuratOBJ, label = TRUE) + NoLegend()

## Vplots for the features selected

DefaultAssay(SeuratOBJ) <- "RNA"

seurat_name <- str_extract(Seurat_base_name, regex("C\\.\\w+"))
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

message("Plots Completed!")
