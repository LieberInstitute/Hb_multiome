########################################################################
## Rename Seurat with final cell types 
## Save RDS with Seurat renamed idents and table with merged clusters
##
## Authors. CSC
## Date. Jun 30, 2025
##
## Recommended resources on interactive mode: srun --pty --mem=80GB --x11 bash
########################################################################

library("Seurat")
library("dplyr")
library("stringr")
library("here")

# directories

Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
inputSeuratRDS <- here(
    "processed-data", 
    "05_Clustering_ARCr", 
    "05_rename_idents", 
    Seurat_base_name)

input_ct_summary_CSV <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "16_rename_idents_post_mean_ratio_HClust",
    "WNN_full_annotation_meta_data.csv")

processedDir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "17_wnn_hierarchical_clustering_final_ct"
)

## Check directories
if (!dir.exists(processedDir)) {
    dir.create(processedDir)
}

#===============================================================================

message(Sys.time(), " - Load Seurat with WNN clusters")

SeuratOBJ <- readRDS(inputSeuratRDS)
SeuratOBJ
DefaultAssay(SeuratOBJ) <- "RNA"
class(SeuratOBJ[["ATAC"]])
current_ct <- levels(SeuratOBJ)
current_ct <- sort(current_ct)

#===============================================================================


message(Sys.time(), " - Load full annotation to renmmae Seurat idents with final cell-types")

summary_ct_df <- read.csv(input_ct_summary_CSV)
colnames(summary_ct_df)

## verify clustering order and ct assignation
summary_ct_df_sorted <- summary_ct_df |>
    arrange(cluster) |>
    select(cluster, ct_final)
summary_ct_df_sorted
## remove (*) and compose and update ct    
new_ct <- paste0(summary_ct_df$cluster, ".", summary_ct_df$ct_final)
new_ct <- gsub("\\*", "", new_ct)

## verify
cell_types_df <- data.frame(
    old_ct = current_ct,
    new_ct = new_ct
)
cell_types_df
# > cell_types_df
#           old_ct          new_ct
# 1    C.01.Inhib.Thal C.01.Inhib.Thal
# 2         C.02.Oligo      C.02.Oligo
# 3    C.03.Excit.Thal C.03.Excit.Thal
# 4    C.04.Excit.Thal      C.04.LHb.4
# 5       C.05.LHb.2.7    C.05.LHb.2.7
# 6  C.06.ExcitT.LHb.4      C.06.LHb.4
# 7         C.07.MHb.2      C.07.MHb.2
# 8         C.08.LHb.4      C.08.LHb.4
# 9  C.09.ExcitT.LHb.4      C.09.LHb.4


names(new_ct) <- current_ct
## rename idents with new cell-types
SeuratOBJ <- RenameIdents(SeuratOBJ, new_ct)
levels(SeuratOBJ)
# [1] "C.01.Inhib.Thal" "C.02.Oligo"      "C.03.Excit.Thal" "C.04.LHb.4"     
# [5] "C.05.LHb.2.7"    "C.06.LHb.4"      "C.07.MHb.2"      "C.08.LHb.4"     
# [9] "C.09.LHb.4"      "C.10.MHb.1"      "C.11.MHb.1.2"    "C.12.Excit.Thal"
# [13] "C.13.LHb.4"      "C.14.MHb.1"      "C.15.Excit.Thal" "C.16.MHb.1.2"   
# [17] "C.17.Excit.Thal" "C.18.LHb.1.3.4"  "C.19.Inhib.Thal" "C.20.Astrocyte" 
# [21] "C.21.Astrocyte"  "C.22.Oligo"      "C.23.LHb.1"      "C.24.LHb.4"     
# [25] "C.25.Excit.Thal" "C.26.OPC"        "C.27.Microglia"  "C.28.Inhib.Thal"
# [29] "C.29.Endo"       "C.30.LHb.7"      "C.31.LHb.4"      "C.32.Excit.Thal"
# [33] "C.33.LHb.1.3"    "C.34.Oligo"      "C.35.Excit.Thal" "C.36.MHb.3"     
# [37] "C.37.Thal"       "C.38.Inhib.Thal" "C.39.Inhib.Thal" "C.40.LHb.4"     
# [41] "C.41.Microglia" 
#head(Idents(SeuratOBJ))

## Set factor levels for identities to arrange clusters, first we want Hb clusters
all_clusters <- as.vector(new_ct)
hb_clusters <- grep("MHb|LHb", all_clusters, value = TRUE)
hb_clusters
# [1] "C.04.LHb.4"     "C.05.LHb.2.7"   "C.06.LHb.4"     "C.07.MHb.2"    
# [5] "C.08.LHb.4"     "C.09.LHb.4"     "C.10.MHb.1"     "C.11.MHb.1.2"  
# [9] "C.13.LHb.4"     "C.14.MHb.1"     "C.16.MHb.1.2"   "C.18.LHb.1.3.4"
# [13] "C.23.LHb.1"     "C.24.LHb.4"     "C.30.LHb.4"     "C.31.LHb.4"    
# [17] "C.33.LHb.1.3"   "C.36.MHb.3"     "C.40.LHb.4"  
no_hb_clust <- grep("MHb|LHb", all_clusters, value = TRUE, invert = TRUE)
no_hb_clust

# ensure all clusters are included
new_levels <- c(hb_clusters, setdiff(all_clusters, hb_clusters))
new_levels
# Apply the new order to Seurat object identities
SeuratOBJ <- SetIdent(SeuratOBJ, value = factor(Idents(SeuratOBJ), levels = new_levels))

message("Added new cluster arrangement")

levels(SeuratOBJ)
unique(Idents(SeuratOBJ))

## =============================================================================

## Remove clusters with wear cell-types

# Drop cluster "C.34.Oligo"
# For details go to: https://github.com/LieberInstitute/Hb_multiome/blob/da5cd6c9ab6c51a63553fbb2cc1080260e089628/processed-data/05_Clustering_ARCr/16_rename_idents_post_mean_ratio_HClust/WNN_full_annotation_meta_data.csv 
cells_to_keep <- WhichCells(SeuratOBJ, idents = NULL)[Idents(SeuratOBJ) != "C.34.Oligo"]
SeuratOBJ <- subset(SeuratOBJ, cells = cells_to_keep)
levels(SeuratOBJ)
unique(Idents(SeuratOBJ))
head(Idents(SeuratOBJ))


## =============================================================================
## Add x meta-cluster as column; eg: MHb, LHb and No-Habenula

## Update "cluster_ann" column to Seurat meta-data for visualizations
colnames(SeuratOBJ@meta.data)

SeuratOBJ$cluster_ann <- Idents(SeuratOBJ)

assign_merged_clusters <- function(cluster_vector, cluster_groups, default = "Other") {
    # cluster_vector : character vector of cluster IDs (e.g. from Idents(SeuratOBJ))
    # cluster_groups : named list, names = merged group names, values = cluster IDs to merge
    # default : fallback category
    
    # Initialize with default value
    merged_vector <- rep(default, length(cluster_vector))
    
    # For each merged group, overwrite matching values
    for (group_name in names(cluster_groups)) {
        matched_idx <- cluster_vector %in% cluster_groups[[group_name]]
        merged_vector[matched_idx] <- group_name
    }
    
    return(merged_vector)
}

## Define new classes to merge

all_Thal_only <- grep("Thal", no_hb_clust, value = TRUE)

cluster_merged_groups <- list(
    LHb = grep("LHb", hb_clusters, value = TRUE),
    MHb = grep("MHb", hb_clusters, value = TRUE)
) |>
    append(list(
        Oligo = grep("Oligo", no_hb_clust, value = TRUE),
        Astrocyte = grep("Astrocyte", no_hb_clust, value = TRUE),
        OPC = grep("OPC", no_hb_clust, value = TRUE),
        Microglia = grep("Microglia", no_hb_clust, value = TRUE),
        Endo = grep("Endo", no_hb_clust, value = TRUE),
        Inhib_Thal = grep("Inhib.Thal", no_hb_clust, value = TRUE),
        Excit_Thal = grep("Excit.Thal", no_hb_clust, value = TRUE),
        Thal = if (length(all_Thal_only) > 0) all_Thal_only[!grepl("\\.Excit\\.Thal|\\.Inhib\\.Thal", all_Thal_only)] else character(0)
    ))

# Add meta-data "merged_cluster" with x merged clusters defined above
current_idents <- as.character(Idents(SeuratOBJ))
merged_cluster <- assign_merged_clusters(current_idents, cluster_merged_groups)
SeuratOBJ$merged_cluster <- merged_cluster

# check
table(SeuratOBJ$merged_cluster)
# Astrocyte       Endo Excit_Thal Inhib_Thal        LHb        MHb  Microglia 
# 2684        343       9738       6024      19673      10944        663 
# Oligo        OPC       Thal 
# 4882        638        111 

## Table summary of merged clusters
merged_table <- table(SeuratOBJ$merged_cluster) |> as.data.frame()
colnames(merged_table) <- c("merged_cluster", "cell_count")
merged_table <- merged_table |> 
    mutate(
        percent = round(100 * cell_count / sum(cell_count), 2)
    )

merged_table <- merged_table |> arrange(desc(cell_count))
print(merged_table)

f_file <- here(processedDir, "merged_clusters_summary.csv")
write.csv(merged_table, file = f_file, row.names = FALSE)

message("Summary table saved at: ", f_file)

message("Clusters rearrenged in big categories merged!")


## =============================================================================

## save RDS
rds_file_name <- here(processedDir, paste0(Seurat_base_name, "_renamed_visium_HD.rds"))
saveRDS(SeuratOBJ, rds_file_name)

message("Seurat with clusters renamed saved!")


# library("slurmjobs")
# job_single(
#     "17_wnn_clustering_final_ct", 
#     cores = 2, 
#     partition = "katun", 
#     memory = "80G", 
#     create_shell = TRUE
#     )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



