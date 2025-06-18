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
inputCVS_cell_types_summary <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze_v3", 
                         "FULL_SUMMARY_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_allTypes_v3.csv")

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


## load summary with cell_types and percentages by clusters 
summary_ct_df <- read.csv(inputCVS_cell_types_summary)
head(summary_ct_df)
colnames(summary_ct_df)
summary_ct_df |>
    select(seurat_clusters, Perc.Cluster, cell_types) |>
    head()
# seurat_clusters Perc.Cluster
# 1               1        7.012
# 2               2        6.052
# 3               3        5.872
# 4               4        5.226
# 5               5        4.975
# 6               6        4.788
# cell_types
# 1 DD_Inhib.Thal (10), DD_LHb (1), LB_excitatory_neuron (1), LB_inhibitory_neuron (2), LB_Thalamus/MDm (1)
# 2                                                     DD_Oligo (5), LB_neuron (2), LB_oligodendrocyte (1)
# 3                                                        DD_Astrocyte (1), DD_Endo (4), DD_Excit.Thal (1)
# 4                DD_Excit.Thal (6), DD_Inhib.Thal (1), DD_LHb (2), DD_OPC (1), LB_LHB neuron specific (1)
# 5                   DD_Excit.Thal (1), DD_LHb (10), LB_Hb neuron specific (2), LB_LHB neuron specific (1)
# 6                                                                           DD_Excit.Thal (6), DD_OPC (1)

tail(summary_ct_df)

# only retain cluster IDs, remove not numeric rows from the summary (total)
summary_ct_df <- summary_ct_df %>%
    filter(grepl("^\\d+$", seurat_clusters))
summary_ct_df$seurat_clusters
# [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "14" "15" "16"
# [16] "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "31" "33"
# [31] "34" "35" "36" "37" "41" "32" "39" "30" "38" "40" "42" "13"

# set pad of 2 digits to clusters
summary_ct_df$cluster <- paste0("C.", sprintf('%02d', as.numeric(as.character(summary_ct_df$seurat_clusters))))
head(summary_ct_df)

# rename and format some columns
summary_ct_df <- summary_ct_df |>
    rename(cluster_percentage = Perc.Cluster, 
           all_cell_types_match = cell_types)  |>
    mutate(cluster_percentage = paste0(cluster_percentage, "%"))
head(summary_ct_df)

## Polished annotation
# Manually curated annotations to annotate Seurat clusters
# Details on: https://github.com/LieberInstitute/Hb_multiome/tree/8b614669db8f7b833e5a793fa01326f839c57cd8/data 
cell_types_curated <- data.frame(
    cell_type = c(
        cell_type_pseudo <- c(
            "undetermined", "DD_Oligo", "undetermined", "undetermined", "DD_LHb",
            "DD_Excit.Thal", "DD_MHb", "undetermined", "undetermined", "DD_MHb",
            "DD_MHb", "undetermined", "no-match", "DD_MHb", "DD_Excit.Thal",
            "DD_MHb", "DD_Excit.Thal", "DD_LHb", "DD_Inhib.Thal", "DD_Astrocyte",
            "DD_Astrocyte", "undetermined", "DD_LHb", "DD_LHb", "undetermined",
            "DD_OPC", "DD_Microglia", "DD_Inhib.Thal", "DD_Endo", "DD_LHb",
            "DD_Excit.Thal", "undetermined", "DD_LHb", "DD_Oligo", "undetermined",
            "DD_MHb", "undetermined", "DD_Inhib.Thal", "DD_Inhib.Thal", "DD_LHb",
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
full_annotation_df <- left_join(cell_types_curated, summary_ct_df, by = "cluster")
colnames(full_annotation_df)

# filtr and sort columns
full_annotation_df <- full_annotation_df[, c("cluster", "cell_type", "cluster_percentage", "all_cell_types_match")]

# ==============================================================================
# save full summary with cluster annotations curated

full_annotation_df
#       cluster      cell_type cluster_percentage
# 1    C.01 undetermined              7.01%
# 2    C.02       DD_Oligo              6.05%
# 3    C.03 undetermined              5.87%
# 4    C.04 undetermined              5.23%
# 5    C.05         DD_LHb              4.97%
# 6    C.06   DD_Excit.Thal              4.79%
# all_cell_types_match
# 1 DD_Inhib.Thal, LB_Thalamus/MDm, LB_inhibitory_neuron, undetermined, DD_LHb, LB_excitatory_neuron
# 2                                            DD_Oligo, undetermined, LB_oligodendrocyte, LB_neuron
# 3                                               undetermined, DD_Endo, DD_Astrocyte, DD_Excit.Thal
# 4               DD_Excit.Thal, undetermined, DD_Inhib.Thal, DD_LHb, DD_OPC, LB_LHB neuron specific
# 5               DD_LHb, LB_LHB neuron specific, undetermined, LB_Hb neuron specific, DD_Excit.Thal
# 6                                                              undetermined, DD_Excit.Thal, DD_OPC


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
head(levels(SeuratOBJ))
head(Idents(SeuratOBJ))

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
# [1] "C.05.DD_LHb"        "C.07.DD_MHb"        "C.10.DD_MHb"       
# [4] "C.11.DD_MHb"        "C.14.DD_MHb"        "C.16.DD_MHb"       
# [7] "C.18.DD_LHb"        "C.23.DD_LHb"        "C.24.DD_LHb"       
# [10] "C.30.DD_LHb"        "C.33.DD_LHb"        "C.36.DD_MHb"       
# [13] "C.40.DD_LHb"        "C.01.undetermined"  "C.02.DD_Oligo"     
# [16] "C.03.undetermined"  "C.04.undetermined"  "C.06.DD_Excit.Thal"
# [19] "C.08.undetermined"  "C.09.undetermined"  "C.12.undetermined" 
# [22] "C.13.no-match"      "C.15.DD_Excit.Thal" "C.17.DD_Excit.Thal"
# [25] "C.19.DD_Inhib.Thal" "C.20.DD_Astrocyte"  "C.21.DD_Astrocyte" 
# [28] "C.22.undetermined"  "C.25.undetermined"  "C.26.DD_OPC"       
# [31] "C.27.DD_Microglia"  "C.28.DD_Inhib.Thal" "C.29.DD_Endo"      
# [34] "C.31.DD_Excit.Thal" "C.32.undetermined"  "C.34.DD_Oligo"     
# [37] "C.35.undetermined"  "C.37.undetermined"  "C.38.DD_Inhib.Thal"
# [40] "C.39.DD_Inhib.Thal" "C.41.DD_Microglia"  "C.42.no-match"  

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

