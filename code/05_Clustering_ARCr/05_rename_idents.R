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
plotDir <- here("plots", "05_Clustering_ARCr", "05_rename_idents")


## read input arguments ( name of RDS Seurat file with wnn clustering to parse )
# Seurat_base_name <- commandArgs(trailingOnly = TRUE)

## For testing: 
Seurat_base_name <- "seurat.norm_counts_Harmony_ARCr_QCed_WNN_k30_C.louvain_lsi_r1"

seurat_RDSname <- here(inputRDS_Dir, paste0(Seurat_base_name, ".rds"))
     
SeuratOBJ <- readRDS(seurat_RDSname)
## verification
length(Cells(x = SeuratOBJ))
nrow(unique(SeuratOBJ[["seurat_clusters"]]))



## (1) Load DEG annotated

file_name <- here(inputCVS_Dir_Ann, "ARCr_QCed_WNN_k30_C.louvain_lsi_r1_DEG_Top50_WNN_DD_LB_matching_markers_integrated.csv")
file_ann <- read.csv(file_name)
head(file_ann)


## Prepare a list with cell-types to rename idents (takes as input the annotated cell-types)

## get unique annotated cell-types (clusters) containing the 'Hb' word
get_hb_words <- as.data.frame(file_ann) |> filter(grepl('Hb', cell_type))
get_hb_words <- unique(get_hb_words[c("cluster", "cell_type")])
## collapse the redundancy clusters
collapsed_cell_types <- aggregate(cell_type ~ cluster, data = get_hb_words, FUN = function(x) paste(x, collapse = ", "))
## set pad of 2 digits to clusters
collapsed_cell_types$cluster <- paste0("C.", sprintf('%02d', collapsed_cell_types$cluster))
v_hb_short_cell_types <- substring(collapsed_cell_types$cell_type, regexpr(" ", collapsed_cell_types$cell_type) + 1)
v_hb_short_cell_types <- substring(v_hb_short_cell_types, 1, 10)
collapsed_cell_types$cell_type <- c(v_hb_short_cell_types)

head(collapsed_cell_types)
# cluster  cell_type
# 1    C.00     DD_MHb
# 2    C.01     DD_LHb
# 3    C.02     DD_MHb
# 4    C.03     DD_LHb
# 5    C.06 LB_Hb neur
# 6    C.08 LB_Hb neur



## (2) extract unique clusters in ascending order and prepare a named list

Seurat_clusterIDS <- as.integer(levels(SeuratOBJ$seurat_clusters))
# Seurat_clusterIDS <- paste0("C.", sprintf('%02d',Seurat_clusterIDS))
# Seurat_clusterIDS <- as.data.frame(Seurat_clusterIDS)
# colnames(Seurat_clusterIDS) = "cluster"

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

v_new_clusterIDS <- paste(Seurat_clusterIDS_new$cluster, Seurat_clusterIDS_new$cell_type)

levels(SeuratOBJ)

SeuratOBJ_renamed <- RenameIdents(object = c(levels(SeuratOBJ)), c(v_new_clusterIDS))

levels(SeuratOBJ)

