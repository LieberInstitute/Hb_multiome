########################################################################
## Search Peaks on ATAC modality: Uses Signac::CallPeaks() to redo ATAC with local peaks (wnn clusters)
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
library("tidyverse")
library("dbplyr")
library("here")


## read input arguments
args = commandArgs(trailingOnly = TRUE)
resolution_level <- args[2]
# resolution_level = "Fine" # 42 clusters
# resolution_level = "Broad"  # 8 cell-types
# resolution_level = "Mid" # 8 cell-types

if (length(resolution_level)) {
    message(
        "Processing peaks for resolution:\n",
        resolution_level
    )
    f_sufix <- paste0(".resolution.", resolution_level)
    f_sufix
} else {
    message("Input arguments missed")
    stop()
}


# Check/create directories
## clusters renamed for Spatial-Registration on Visium project

inputRDS_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "22_add_mid_level_clustering"
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

message("Loading Seurat ... ")

# Use Seurat with final ct-annotations
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))

##==============================================================================
## Verification

# colnames(SeuratOBJ@meta.data)
msg <- switch(resolution_level,
        "Fine" = paste0("Fine-level Cell-Types:\n ", paste(sort(unique(SeuratOBJ$cluster_ann)), collapse = "\n")),
        "Broad" = paste0("Broad-level Cell-Types:\n", paste(sort(unique(SeuratOBJ$merged_cluster)), collapse = "\n")),
        "Mid" = paste0("Mid-level Cell-Types:\n", paste(sort(unique(SeuratOBJ$mid_cluster)), collapse = "\n"))
)

message(msg)
message("Processing ", length(Cells(SeuratOBJ)), " cells")

## set col name / should exists as meta-data
## testing:
# resolution_level="Broad"
# resolution_level="Fine"
meta_col <- case_when(
    resolution_level=="Fine" ~ "cluster_ann",
    resolution_level=="Broad" ~ "merged_cluster",
    resolution_level=="Mid" ~ "mid_cluster"
)


##==============================================================================
## Set desired meta-data as current level

# Set Seurat identities from a metadata column
set_idents_from_meta <- function(seurat_obj, meta_col, level_order = NULL, na_fill = "Unknown") {
    ## double check level exist on meta-data
    if (!meta_col %in% colnames(seurat_obj@meta.data)) {
        stop("Meta column '", meta_col, "' not found in SeuratOBJ@meta.data")
    }
    # extract target vector
    target_vec <- as.character(seurat_obj[[meta_col]][, 1])
    
    # decide levels and keep appearance order
    if (is.null(level_order)) {
        level_order <- sort(unique(target_vec))
    }
    
    # only update if different from current Idents
    current_idents <- as.character(Idents(seurat_obj))
    if (!identical(current_idents, target_vec)) {
        seurat_obj <- SetIdent(seurat_obj, value = factor(target_vec, levels = level_order))
        message("Idents set from meta column '", meta_col, "'.")
    } else {
        message("Idents already match '", meta_col, "', nothing to do.")
    }
    
    return(seurat_obj)
}

SeuratOBJ <- set_idents_from_meta(SeuratOBJ, meta_col = meta_col)

# specific order:
# desired_levels <- c("LHb.1","LHb.1.3","LHb.1.3.4","LHb.2.7","LHb.4","LHb.7",
#                     "MHb.1","MHb.1.2","MHb.2","MHb.3",
#                     "Excit.Thal","Inhib.Thal","Astrocyte","Oligo","OPC","Microglia","Endo","Thal")
# SeuratOBJ <- set_idents_from_meta(SeuratOBJ, "mid_cluster", level_order = desired_levels)

## sanity check
levels(SeuratOBJ)           # new levels
table(Idents(SeuratOBJ))    # counts by new identities


##==============================================================================

message("Preparing for CallPeaks ...")

DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])
SeuratOBJ[["ATAC"]]
# ChromatinAssay data with 262951 features for 55702 cells
# Variable features: 249866 
# Genome: 
#     Annotation present: TRUE 
# Motifs present: FALSE 
# Fragment files: 10 


##==============================================================================
# call peaks on a single-cell ATAC-seq dataset using MACS2
# peaks will be called independently on each group of cells and then combined

message("Calling peaks at FINE level ... ")

##==============================================================================
# # testing with small cluster: Subset Seurat object where cluster name = C.41.Microglia
# unique(Idents(SeuratOBJ))
# Seurat_subset <- subset(SeuratOBJ, idents = "C.41.Microglia")
# length(Cells(Seurat_subset))
##==============================================================================

## call macs
peaks <- CallPeaks(
    object = SeuratOBJ,
    #object = Seurat_subset,
    group.by = meta_col,
    macs2.path = "/users/csoto/.conda/envs/macs2_conda3_env/bin/macs2",
    verbose = TRUE
)

message("Calling peaks done!")

head(peaks)

# Convert GRanges to data frame and save for further analysis
df_peaks <- as.data.frame(peaks)
head(df_peaks)

f_name <- paste0("macs_peaks_", resolution_level, "_resolution.csv")
write.csv(
    df_peaks, 
    here(outputCSV_Dir, f_name), 
    row.names = FALSE
)

message("Peaks file saved!")

message("All done!!!")


# library("slurmjobs")
# job_single(
#   "01_call_peaks_MACS2",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "80G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript -e \"options(width = 120); sessioninfo::session_info()\"",
#   create_logdir = TRUE
# )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()




# ##==============================================================================
# 
# # Get current cluster identities
# current_idents <- as.character(Idents(seurat_obj))
# 
# # Replace LHb cluster names with 'LHb_merged'
# merged_idents <- ifelse(current_idents %in% LHb_clusters_to_merge, 
#                         "LHb_merged", 
#                         current_idents)
# 
# # Assign new identities to the Seurat object
# Idents(seurat_obj) <- merged_idents
# 
# ##==============================================================================
# 
# 
# # Cell Ranger peaks
# DefaultAssay(SeuratOBJ) <- "ATAC"
# p1 <- CoveragePlot(
#     object = SeuratOBJ,
#     region = gene,
#     features = gene,
#     extend.upstream = up,
#     extend.downstream = down,
#     peaks = TRUE,
#     links = FALSE,
#     annotation = TRUE
# ) + ggtitle("Cell Ranger peaks")
# 
# # MACS2 peaks
# DefaultAssay(SeuratOBJ) <- "ATAC_MACS2"
# p2 <- CoveragePlot(
#     object = SeuratOBJ,
#     region = gene,
#     features = gene,
#     extend.upstream = up,
#     extend.downstream = down,
#     peaks = TRUE,
#     links = FALSE,
#     annotation = TRUE
# ) + ggtitle("MACS2 peaks (pseudobulk)")
# 
# # Combine
# p_combined <- p1 / p2

# ggsave(
#     filename = file.path(plotDir, paste0("coverage_before_after_MACS2_", gene, ".pdf")),
#     plot = p_combined, width = 10, height = 6
# )



