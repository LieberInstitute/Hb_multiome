########################################################################
## Annotate WNN Clusters: all cluster with cell-types identified with the human pilot annotation
## INPUT:
##      (1) GEX DEG annotation
##      (2) Seurat with WNN 
## OUPUT:
##      (1) Seurats with annotated idents and MERGED meta-data added
##      (2) Full summary with cluster annotations curated (supplementary material)
##      (3) Some visualizations: VPlots, DimPlot, Feature, DotPlot ...
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
outputCSV_Dir <- here("data")
plotDir <- here("plots", "05_Clustering_ARCr", "05_rename_idents")

## Check directories
if (!dir.exists(plotDir)) {dir.create(plotDir)}
if (!dir.exists(outputRDS_Dir)) {dir.create(outputRDS_Dir)}

## read input arguments ( name of RDS Seurat file with wnn clustering to parse )
Seurat_base_name <- commandArgs(trailingOnly = TRUE)

Seurat_base_name <- Seurat_base_name[[2]]
## For debugging: 
# Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2"

message("Processing ", Seurat_base_name)


##### (1) Load Seurat with WNN idents given by default 

seurat_RDSname <- here(inputRDS_Dir, paste0(Seurat_base_name, ".rds"))
     
SeuratOBJ <- readRDS(seurat_RDSname)

total_cells <- length(Cells(x = SeuratOBJ))

message("Renaming ", nrow(unique(SeuratOBJ[["seurat_clusters"]])), " clusters for ", Seurat_base_name)

##### (2) Load DEG with ident annotation

## extract a shorter name to save files 
tmp_wd <- str_extract(Seurat_base_name, pattern = "WNN\\w*\\.\\w*")
deg_file <- paste0(tmp_wd, "_cellTypes_integrated_top50.csv")
# WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r1_cellTypes_integrated_top50.csv
file_name <- here(inputCVS_Dir_Ann, deg_file)
file_ann <- read.csv(file_name)
head(file_ann)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster       gene cell_type
# 1     0  -3.230185 0.222 0.693         0       0      CPNE4      <NA>
# 2     0  -3.687100 0.184 0.640         0       0       VAV3      <NA>
# 3     0  -4.194794 0.029 0.460         0       0 AC119673.2      <NA>
length(unique(Idents(SeuratOBJ))) == length(unique(file_ann$cluster))


## Prepare a list with cell-types to rename idents (takes as input the annotated cell-types)

## get unique annotated cell-types (clusters) containing the 'Hb' word
# cell_types_df <- as.data.frame(file_ann) |> filter(grepl('Hb', cell_type))
cell_types_df <- as.data.frame(file_ann)

## replace NA values with the string "undeterminated" 
cell_types_df <- cell_types_df |>
    mutate(cell_type_longer = ifelse(is.na(cell_type), "undeterminated", cell_type))
table(is.na(cell_types_df$cell_type_longer))

cell_types_df <- unique(cell_types_df[c("cluster", "cell_type_longer")])
head(cell_types_df)
unique(cell_types_df$cluster)

## Calculate % of cells by cluster

f_get_percent_label <- function(cl) {
  ## for testing:  cl = 0
  clust_size <- length(SeuratOBJ$seurat_clusters[SeuratOBJ[["seurat_clusters"]]==cl])
  ## longer label
  # clust_label <- paste0("[", clustSize = clust_size, " / ", paste0(round(((clust_size*100) / total_cells), 2), "%"), "]")
  ## shorter label
  clust_label <- paste0(round(((clust_size*100) / total_cells), 2), "%")
  return(clust_label)               
}

##  add labels to df: cluster id, short label, long_label, percentages, etc

v_unique_hb_clust <- unlist(unique(cell_types_df["cluster"]))
cells_percents <- map(v_unique_hb_clust, ~ f_get_percent_label(.x))
head(cells_percents)

cell_types_with_percents_df <- data.frame(
  cluster = c(v_unique_hb_clust),
  percent = c(unlist(cells_percents)))
cell_types_with_percents_df
#     cluster percent
# cluster1        1   7.01%
# cluster2        2   6.05%
# cluster3        3   5.87%

# collapse the redundancy clusters
collapsed_cell_types <- aggregate(cell_type_longer ~ cluster, data = cell_types_df, FUN = function(x) paste(x, collapse = ", "))

# set pad of 2 digits to clusters
collapsed_cell_types$cluster <- paste0("C.", sprintf('%02d', collapsed_cell_types$cluster))

# shorter the cell_type label
# pseudo annotation using the first match == ONLY for EDA
collapsed_cell_types$cell_type_pseudo <- sub(",.*", "_pseudo", collapsed_cell_types$cell_type)
collapsed_cell_types$cluster_percentage <- cell_types_with_percents_df$percent

## Polished annotation
# Manually curated annotations to annotate Seurat clusters
# Details on: https://github.com/LieberInstitute/Hb_multiome/tree/8b614669db8f7b833e5a793fa01326f839c57cd8/data 
cell_types_curated <- data.frame(
    cell_type = c(
        cell_type_pseudo <- c(
            "undeterminated", "DD_Oligo", "undeterminated", "undeterminated", "DD_LHb",
            "DD_Exit.Thal", "DD_MHb", "undeterminated", "undeterminated", "DD_MHb",
            "DD_MHb", "undeterminated", "no-match", "DD_MHb", "DD_Exit.Thal",
            "DD_MHb", "DD_Exit.Thal", "DD_LHb", "DD_Inhib.Thal", "DD_Astrocyte",
            "DD_Astrocyte", "undeterminated", "DD_LHb", "DD_LHb", "undeterminated",
            "DD_OPC", "DD_Microglia", "DD_Inhib.Thal", "DD_Endo", "DD_LHb",
            "DD_Exit.Thal", "undeterminated", "DD_LHb", "DD_Oligo", "undeterminated",
            "DD_MHb", "undeterminated", "DD_Inhib.Thal", "DD_Inhib.Thal", "DD_LHb",
            "DD_Microglia", "no-match"
        )
    ),
    cluster = c(
        "C.01", "C.02", "C.03", "C.04", "C.05", "C.06", "C.07", "C.08", "C.09", "C.10",
        "C.11", "C.12", "C.13", "C.14", "C.15", "C.16", "C.17", "C.18", "C.19", "C.20",
        "C.21", "C.22", "C.23", "C.24", "C.25", "C.26", "C.27", "C.28", "C.29", "C.30",
        "C.31", "C.32", "C.33", "C.34", "C.35", "C.36", "C.37", "C.38", "C.39", "C.40",
        "C.41", "C.42"
    ),
    stringsAsFactors = FALSE  # Optional, prevents conversion to factors
)
tail(cell_types_curated)

## inner join pseudo annotation (EDA) with curated annotation for supplemental material
full_annotation_df <- left_join(cell_types_curated, collapsed_cell_types, by = "cluster")
colnames(full_annotation_df)
# order columns
colnames(full_annotation_df)[colnames(full_annotation_df) == "cell_type_longer"] <- "all_cell_types_match"
full_annotation_df <- full_annotation_df[, c("cluster", "cell_type", "cluster_percentage", "all_cell_types_match", "cell_type_pseudo")]


# ==============================================================================
# save full summary with cluster annotations curated

full_annotation_df
#       cluster      cell_type cluster_percentage
# 1    C.01 undeterminated              7.01%
# 2    C.02       DD_Oligo              6.05%
# 3    C.03 undeterminated              5.87%
# 4    C.04 undeterminated              5.23%
# 5    C.05         DD_LHb              4.97%
# 6    C.06   DD_Exit.Thal              4.79%
# all_cell_types_match
# 1 DD_Inhib.Thal, LB_Thalamus/MDm, LB_inhibitory_neuron, undeterminated, DD_LHb, LB_excitatory_neuron
# 2                                            DD_Oligo, undeterminated, LB_oligodendrocyte, LB_neuron
# 3                                               undeterminated, DD_Endo, DD_Astrocyte, DD_Excit.Thal
# 4               DD_Excit.Thal, undeterminated, DD_Inhib.Thal, DD_LHb, DD_OPC, LB_LHB neuron specific
# 5               DD_LHb, LB_LHB neuron specific, undeterminated, LB_Hb neuron specific, DD_Excit.Thal
# 6                                                              undeterminated, DD_Excit.Thal, DD_OPC


write.csv(full_annotation_df, here(outputCSV_Dir, "full_annotation_meta_data.csv"), row.names = FALSE)
# ==============================================================================


## extract original idents (clusters) and prepare the new ident names to rename Seurat clusters

# sanity check if any Seurat clusters was excluded
Seurat_clusterIDS <- as.integer(levels(SeuratOBJ$seurat_clusters))
Seurat_clusterIDS <- paste0("C.", sprintf('%02d',Seurat_clusterIDS))
if (!identical(Seurat_clusterIDS, full_annotation_df$cluster)) {
    stop("Cluster IDs in Seurat object do not match those in the annotation data frame.")
}

# Extract numeric part from "C.XX" and convert to character
full_annotation_df$cluster_id <- as.character(as.numeric(sub("C\\.", "", full_annotation_df$cluster)))

# new_names = current cluster IDs in SeuratOBJ + full_annotation_df$cell_type (CURATED ANNOTATION)
new_names <- setNames(paste0(full_annotation_df$cluster, ".", full_annotation_df$cell_type), full_annotation_df$cluster_id)
new_names

## rename idents 
# Idents(SeuratOBJ)
SeuratOBJ <- RenameIdents(object = SeuratOBJ, new_names)

# verification
levels(SeuratOBJ)
# [1] "C.01.undeterminated" "C.02.DD_Oligo"       "C.03.undeterminated"
# [4] "C.04.undeterminated" "C.05.DD_LHb"         "C.06.DD_Exit.Thal"  
# [7] "C.07.DD_MHb"         "C.08.undeterminated" "C.09.undeterminated"
# [10] "C.10.DD_MHb"         "C.11.DD_MHb"         "C.12.undeterminated"
# [13] "C.13.no-match"       "C.14.DD_MHb"         "C.15.DD_Exit.Thal"  
# [16] "C.16.DD_MHb"         "C.17.DD_Exit.Thal"   "C.18.DD_LHb"        
# [19] "C.19.DD_Inhib.Thal"  "C.20.DD_Astrocyte"   "C.21.DD_Astrocyte"  
# [22] "C.22.undeterminated" "C.23.DD_LHb"         "C.24.DD_LHb"        
# [25] "C.25.undeterminated" "C.26.DD_OPC"         "C.27.DD_Microglia"  
# [28] "C.28.DD_Inhib.Thal"  "C.29.DD_Endo"        "C.30.DD_LHb"        
# [31] "C.31.DD_Exit.Thal"   "C.32.undeterminated" "C.33.DD_LHb"        
# [34] "C.34.DD_Oligo"       "C.35.undeterminated" "C.36.DD_MHb"        
# [37] "C.37.undeterminated" "C.38.DD_Inhib.Thal"  "C.39.DD_Inhib.Thal" 
# [40] "C.40.DD_LHb"         "C.41.DD_Microglia"   "C.42.no-match"
head(Idents(SeuratOBJ))
# S04_AAACAGCCAGAATGAC-1 S04_AAACAGCCAGCAAGGC-1 S04_AAACATGCACCTGGTG-1 
# C.25.undeterminated    C.04.undeterminated    C.09.undeterminated 
# S04_AAACATGCAGGATGGC-1 S04_AAACATGCAGTAATAG-1 S04_AAACATGCATAAGTCT-1 
# C.04.undeterminated    C.01.undeterminated    C.01.undeterminated 
# 42 Levels: C.01.undeterminated C.02.DD_Oligo ... C.42.no-match

# Set new factor levels for identities to first plot Hb clusters
all_clusters <- levels(SeuratOBJ)
hb_clusters <- grep("MHb|LHb", all_clusters, value = TRUE)
hb_clusters
# [1] "C.05 DD_LHb" "C.07 DD_MHb" "C.10 DD_MHb" "C.11 DD_MHb" "C.14 DD_MHb"
# [6] "C.16 DD_MHb" "C.18 DD_LHb" "C.23 DD_LHb" "C.24 DD_LHb" "C.30 DD_LHb"
# [11] "C.33 DD_LHb" "C.36 DD_MHb" "C.40 DD_LHb"
no_hb_clust <- grep("MHb|LHb", all_clusters, value = TRUE, invert = TRUE)
no_hb_clust

# ensure all clusters are included
new_levels <- c(hb_clusters, setdiff(all_clusters, hb_clusters))
new_levels
# Apply the new order to Seurat object identities
SeuratOBJ <- SetIdent(SeuratOBJ, value = factor(Idents(SeuratOBJ), levels = new_levels))
levels(SeuratOBJ)
# [1] "C.05.DD_LHb"         "C.07.DD_MHb"         "C.10.DD_MHb"        
# [4] "C.11.DD_MHb"         "C.14.DD_MHb"         "C.16.DD_MHb"        
# [7] "C.18.DD_LHb"         "C.23.DD_LHb"         "C.24.DD_LHb"        
# [10] "C.30.DD_LHb"         "C.33.DD_LHb"         "C.36.DD_MHb"        
# [13] "C.40.DD_LHb"         "C.01.undeterminated" "C.02.DD_Oligo"      
# [16] "C.03.undeterminated" "C.04.undeterminated" "C.06.DD_Exit.Thal"  
# [19] "C.08.undeterminated" "C.09.undeterminated" "C.12.undeterminated"
# [22] "C.13.no-match"       "C.15.DD_Exit.Thal"   "C.17.DD_Exit.Thal"  
# [25] "C.19.DD_Inhib.Thal"  "C.20.DD_Astrocyte"   "C.21.DD_Astrocyte"  
# [28] "C.22.undeterminated" "C.25.undeterminated" "C.26.DD_OPC"        
# [31] "C.27.DD_Microglia"   "C.28.DD_Inhib.Thal"  "C.29.DD_Endo"       
# [34] "C.31.DD_Exit.Thal"   "C.32.undeterminated" "C.34.DD_Oligo"      
# [37] "C.35.undeterminated" "C.37.undeterminated" "C.38.DD_Inhib.Thal" 
# [40] "C.39.DD_Inhib.Thal"  "C.41.DD_Microglia"   "C.42.no-match"   

message("Clusters sorted done!")


## =============================================================================
## Add 3 meta-cluster as column: MHb, LHb and No-Habenula

# extract the ident ID for the 3 meta-groups
LHb_clusters_to_merge <- grep("LHb", hb_clusters, value = TRUE)
MHb_clusters_to_merge <- grep("MHb", hb_clusters, value = TRUE)
LHb_clusters_to_merge
# [1] "C.05.DD_LHb" "C.18.DD_LHb" "C.23.DD_LHb" "C.24.DD_LHb" "C.30.DD_LHb"
# [6] "C.33.DD_LHb" "C.40.DD_LHb"
MHb_clusters_to_merge
# [1] "C.07.DD_MHb" "C.10.DD_MHb" "C.11.DD_MHb" "C.14.DD_MHb" "C.16.DD_MHb"
# [6] "C.36.DD_MHb"
no_hb_clust

# Add meta-data "merged_cluster" with 3 merged clusters classes: LHb, MHb and No-Hb clusters
current_idents <- as.character(Idents(SeuratOBJ))
# Assign merged labels
merged_cluster <- ifelse(current_idents %in% LHb_clusters_to_merge, "LHb_merged",
                         ifelse(current_idents %in% MHb_clusters_to_merge, "MHb_merged",
                                ifelse(current_idents %in% no_hb_clust, "No-Hb_merged", current_idents)))
unique(merged_cluster)

# add to new metadata MERGED ident labels for further analysis
SeuratOBJ$merged_cluster <- merged_cluster
unique(SeuratOBJ$merged_cluster)
#[1] "No-Hb_merged" "MHb_merged"   "LHb_merged" 

## save RDS
rds_file_name <- here(outputRDS_Dir, paste0(Seurat_base_name, "_renamed_visium.rds"))
saveRDS(SeuratOBJ, rds_file_name)

message("New seurat with clusters annotated and `merged_cluster` meta-data saved!")


## =============================================================================
## Some visualizations: VPlots, DimPlot, Feature, DotPlot ...

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
#Reductions(SeuratOBJ)
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

