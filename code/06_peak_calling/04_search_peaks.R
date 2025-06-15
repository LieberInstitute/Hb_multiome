########################################################################
## Search Peaks on ATAC modality
##
## Authors. CSC
## Date. June 13, 2025
## Recommended resources on interactive mode: srun --pty --mem=80GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
## Resources:
## https://stuartlab.org/signac/articles/peak_calling
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
library("here")

# Check/create directories
## clusters renamed for Spatial-Registration on Visium project

inputRDS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "05_rename_idents"
)
plotDir <- here(
    "plots",
    "06_peak_calling"
)
outputCSV_Dir <- here("processed-data", "06_peak_calling")

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(outputCSV_Dir)) {
    dir.create(outputCSV_Dir)
}

## Load Seurat
# Use Seurat with clusters renamed, also used for Spatial-Registration on Visium project

Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
#Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
levels(SeuratOBJ)

DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])
SeuratOBJ[["ATAC"]]
# ChromatinAssay data with 262951 features for 55702 cells
# Variable features: 249866 
# Genome: 
#     Annotation present: TRUE 
# Motifs present: FALSE 
# Fragment files: 10 

# filter Hb clusters
Seurat_subset <- subset(SeuratOBJ, idents = grep("MHb|LHb", levels(SeuratOBJ), value = TRUE))
levels(Seurat_subset)
# [1] "C.05.DD_LHb" "C.07.DD_MHb" "C.10.DD_MHb" "C.11.DD_MHb" "C.14.DD_MHb"
# [6] "C.16.DD_MHb" "C.18.DD_LHb" "C.23.DD_LHb" "C.24.DD_LHb" "C.30.DD_LHb"
# [11] "C.33.DD_LHb" "C.36.DD_MHb" "C.40.DD_LHb"

colnames(Seurat_subset@meta.data)
Seurat_subset@meta.data$seurat_clusters 
Seurat_subset@meta.data$C.leiden_wnn

##==============================================================================
# call peaks on a single-cell ATAC-seq dataset using MACS2

peaks <- CallPeaks(
    object = Seurat_subset,
    group.by = "seurat_clusters"
)

# Convert GRanges to data frame and save as csv file
df_peaks <- as.data.frame(peaks_by_celltype)
head(df_peaks)
write.csv(df_peaks, here(outputCSV_Dir, "Hb_celltype_peaks.csv"), row.names = FALSE)


##==============================================================================

# ## extract base_name to save plots 
# Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+")) 
# Seurat_base_name
# # C.leiden_lsi_r2_renamed_visium
# 
# message("Searching peaks in ATAC modality ...")
# message("Searching on Habenula ATAC clusters ...")
# 
# # Store peak sets per cluster
# cluster_peaks <- list()
# 
# for (clus in clusters) {
#     message("Processing cluster: ", clus)
#     
#     # Subset cells from one cluster
#     cells_in_cluster <- WhichCells(SeuratOBJ, idents = clus)
#     
#     # Subset ATAC counts matrix
#     atac_counts <- GetAssayData(SeuratOBJ, assay = "ATAC", slot = "counts")[, cells_in_cluster]
#     
#     # Find peaks with any signal in this cluster
#     peaks_present <- rownames(atac_counts)[Matrix::rowSums(atac_counts) > 0]
#     
#     # Store
#     cluster_peaks[[clus]] <- peaks_present
# }


# ##==============================================================================
# ## Preparing to compute pair peaks comparison in the Habenula clusters
# 
# all_cluster_IDs <- levels(SeuratOBJ)
# ## extract Hb clusters
# Hb_cluster_IDs <- cluster_IDs[grepl("MHb|LHb", all_cluster_IDs)]
# Hb_cluster_IDs
# # [1] "C.05.DD_LHb" "C.07.DD_MHb" "C.10.DD_MHb" "C.11.DD_MHb" "C.14.DD_MHb"
# # [6] "C.16.DD_MHb" "C.18.DD_LHb" "C.23.DD_LHb" "C.24.DD_LHb" "C.30.DD_LHb"
# # [11] "C.33.DD_LHb" "C.36.DD_MHb" "C.40.DD_LHb"
# pairwise_combinations <- combn(Hb_cluster_IDs, 2, simplify = FALSE)
# pairwise_combinations
# pairwise_df <- do.call(rbind, pairwise_combinations)
# colnames(pairwise_df) <- c("cluster_1", "cluster_2")
# pairwise_df <- as.data.frame(pairwise_df)
# pairwise_df
# 
# 
# for (clust_p in pairwise_df) {
#     da_peaks <- FindMarkers(
#         object = seurat_atac,
#         ident.1 = clust_p["cluster1"],
#         ident.2 = clust_p["cluster2"],
#         test.use = 'LR',
#         min.pct = 0.05
#     )
# }

# library("slurmjobs")
# job_single(
#   "03_coverage_3columns.R",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "60G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript -e \"options(width = 120); sessioninfo::session_info()\"",
#   create_logdir = TRUE
# )


