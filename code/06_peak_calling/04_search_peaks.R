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

# verification
length(Cells(SeuratOBJ))
colnames(SeuratOBJ@meta.data)

all_clusters_df <- data.frame(
    cluster_ann = unique(SeuratOBJ@meta.data[c("cluster_ann")]),
    seurat_cluster = unique(SeuratOBJ@meta.data[c("seurat_clusters")]),
    stringsAsFactors = FALSE  
)
all_clusters_df <- head(all_clusters_df)
print(all_clusters_df, row.names = FALSE)
# cluster_ann           seurat_clusters
# C.25.undetermined              25
# C.04.undetermined               4
# C.09.undetermined               9
# C.01.undetermined               1
# C.13.no-match                   13
# C.08.undetermined               8


##==============================================================================
# call peaks on a single-cell ATAC-seq dataset using MACS2
# I am using pseudo-annotated clusters

peaks <- CallPeaks(
    object = SeuratOBJ,
    group.by = "cluster_ann",
    macs2.path = "/users/csoto/.conda/envs/macs2_conda3_env/bin/macs2"
)
head(peaks)
# GRanges object with 6 ranges and 6 metadata columns:
#     seqnames      ranges strand |            name     score fold_change
# <Rle>   <IRanges>  <Rle> |     <character> <integer>   <numeric>
# [1] GL000194.1 24142-25380      * | Habenula_peak_1        55     1.71652
# [2] GL000194.1 56009-56595      * | Habenula_peak_2        37     1.58937
# [3] GL000194.1 58210-58796      * | Habenula_peak_3        66     1.78009
# [4] GL000194.1 59643-60256      * | Habenula_peak_4        58     1.72924
# [5] GL000194.1 67999-68873      * | Habenula_peak_5        60     1.74195
# [6] GL000194.1 71745-73160      * | Habenula_peak_6        62     1.75467
# neg_log10pvalue_summit neg_log10qvalue_summit relative_summit_position
# <numeric>              <numeric>                <integer>
# [1]                8.62229                5.59984                     1092
# [2]                6.32284                3.73091                      456
# [3]                9.89004                6.62440                      189
# [4]                8.86969                5.80120                      522
# [5]                9.12019                6.00435                      783
# [6]                9.37376                6.20920                      617
# -------
#     seqinfo: 30 sequences from an unspecified genome; no seqlengths

# Convert GRanges to data frame and save as csv file
df_peaks <- as.data.frame(peaks)
head(df_peaks)
write.csv(df_peaks, here(outputCSV_Dir, "Hb_peaks_by_cluster.csv"), row.names = FALSE)


##==============================================================================

# filter peaks for Hb clusters
target_clusters <- grep("MHb|LHb", levels(SeuratOBJ), value = TRUE)
message("Target Hb clusters")
print(target_clusters)
# [1] "C.05.DD_LHb" "C.07.DD_MHb" "C.10.DD_MHb" "C.11.DD_MHb" "C.14.DD_MHb"
# [6] "C.16.DD_MHb" "C.18.DD_LHb" "C.23.DD_LHb" "C.24.DD_LHb" "C.30.DD_LHb"
# [11] "C.33.DD_LHb" "C.36.DD_MHb" "C.40.DD_LHb"

## Subset Peaks where any target cluster is in peak_called_in
matches <- sapply(target_clusters, function(cl) grepl(paste0("\\b", cl, "\\b"), peaks$peak_called_in))
matched_rows <- rowSums(matches) > 0
peaks_subset <- peaks[matched_rows]
peaks_subset <- as.data.frame(peaks_subset)

write.csv(peaks_subset, here(outputCSV_Dir, "peaks_DD_LHb_MHb.csv"), row.names = FALSE)


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
#   "04_search_peaks",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "80G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript -e \"options(width = 120); sessioninfo::session_info()\"",
#   create_logdir = TRUE
# )


