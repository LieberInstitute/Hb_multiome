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
    "17_wnn_clustering_final_ct"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "01_call_peaks_MACS2"
)
outputCSV_Dir <- here(
    "processed-data", 
    "06_peak_calling",
    "01_call_peaks_MACS2"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(outputCSV_Dir)) {
    dir.create(outputCSV_Dir)
}

## Load Seurat
# Use Seurat with final ct-annotations
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
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


message("Processing Peaks for WNN clusters ... ")


all_clusters_df <- data.frame(
    cluster_ann = unique(SeuratOBJ@meta.data[c("cluster_ann")]),
    seurat_cluster = unique(SeuratOBJ@meta.data[c("seurat_clusters")]),
    stringsAsFactors = FALSE  
)
all_clusters_df <- head(all_clusters_df)
print(all_clusters_df, row.names = FALSE)


##==============================================================================
# call peaks on a single-cell ATAC-seq dataset using MACS2
# peaks will be called independently on each group of cells and then combined

message("Calling peaks grouped by wnn cluster ... ")

peaks <- CallPeaks(
    object = SeuratOBJ,
    group.by = "cluster_ann",
    macs2.path = "/users/csoto/.conda/envs/macs2_conda3_env/bin/macs2"
)

message("Calling peaks done!")

head(peaks)

# Convert GRanges to data frame and save for further analysis
df_peaks <- as.data.frame(peaks)
head(df_peaks)
write.csv(
    df_peaks, 
    here(outputCSV_Dir, "all_peaks_by_cluster.csv"), 
    row.names = FALSE
)

# Cell Ranger peaks
DefaultAssay(SeuratOBJ) <- "ATAC"
p1 <- CoveragePlot(
    object = SeuratOBJ,
    region = gene,
    features = gene,
    extend.upstream = up,
    extend.downstream = down,
    peaks = TRUE,
    links = FALSE,
    annotation = TRUE
) + ggtitle("Cell Ranger peaks")

# MACS2 peaks
DefaultAssay(SeuratOBJ) <- "ATAC_MACS2"
p2 <- CoveragePlot(
    object = SeuratOBJ,
    region = gene,
    features = gene,
    extend.upstream = up,
    extend.downstream = down,
    peaks = TRUE,
    links = FALSE,
    annotation = TRUE
) + ggtitle("MACS2 peaks (pseudobulk)")

# Combine
p_combined <- p1 / p2

# ggsave(
#     filename = file.path(plotDir, paste0("coverage_before_after_MACS2_", gene, ".pdf")),
#     plot = p_combined, width = 10, height = 6
# )

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


