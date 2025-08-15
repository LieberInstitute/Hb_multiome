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

message("Loading Seurat ... ")

# Use Seurat with final ct-annotations
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))

# ##==============================================================================
# ## Adding mid level clustering resolution
# 
# # Define mid-level
# clusters_to_merge <- c(SeuratOBJ$cluster_ann)
# 
# # Create a vector mapping old cluster names to MID-LEVEL clusters
# cluster_map <- sapply(clusters_to_merge, function(x) {
#     if (grepl("LHb.4", x)) { return("LHb.4")
#     } else if (grepl("LHb.2.7", x)) { return("LHb.2.7")
#     } else if (grepl("LHb.1$", x)) { return("LHb.1")
#     } else if (grepl("LHb.1.3$", x)) { return("LHb.1.3")
#     } else if (grepl("LHb.1.3.4", x)) { return("LHb.1.3.4")
#     } else if (grepl("LHb.7", x)) { return("LHb.7")
#     } else if (grepl("MHb.1$", x)) { return("MHb.1")
#     } else if (grepl("MHb.1.2", x)) { return("MHb.1.2")
#     } else if (grepl("MHb.2", x)) { return("MHb.2")
#     } else if (grepl("MHb.3", x)) { return("MHb.3")
#     } else if (grepl("Inhib.Thal", x)) { return("Inhib.Thal")
#     } else if (grepl("Excit.Thal", x)) { return("Excit.Thal")
#     } else if (grepl("Astrocyte", x)) { return("Astrocyte")
#     } else if (grepl("OPC", x)) { return("OPC")
#     } else if (grepl("Oligo", x)) { return("Oligo")
#     } else if (grepl("Microglia", x)) { return("Microglia")
#     } else if (grepl("Thal", x)) { return("Thal")
#     } else if (grepl("Endo", x)) { return("Endo")
#     } else { return(NA_character_)  # unmatched case
#     }
# })
# 
# names(cluster_map) <- clusters_to_merge
# # map current identities to LHb/MHb
# mid_idents <- plyr::revalue(as.character(Idents(SeuratOBJ)), cluster_map)
# # assign new identities
# Idents(SeuratOBJ) <- mid_idents
# ## add as meta=data
# SeuratOBJ$mid_cluster <- mid_idents
# # verify mid-levels
# 
# if (length(unique(SeuratOBJ$cluster_ann[is.na(SeuratOBJ$mid_cluster)])) > 0) {
#     stop(
#         "All clusters should be assigned to a mid-level resolution\n",
#         "Levels not assigned\n",
#         paste(unique(SeuratOBJ$cluster_ann[is.na(SeuratOBJ$mid_cluster)]), collapse = "\n")
#     )
# } else {
#     message("Mid-level added!")
#     data.frame(mid_cluster = sort(unique(SeuratOBJ$mid_cluster)))
#     nrow(data.frame(mid_cluster = sort(unique(SeuratOBJ$mid_cluster))))
# }

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
# testing with small cluster: Subset Seurat object where cluster name = C.41.Microglia
unique(Idents(SeuratOBJ))
Seurat_subset <- subset(SeuratOBJ, idents = "C.41.Microglia")
##==============================================================================

length(Cells(Seurat_subset))
#Seurat_subset <- subset(SeuratOBJ, idents = grep("MHb|LHb", Idents(SeuratOBJ), value = TRUE))

peaks <- CallPeaks(
    #object = SeuratOBJ,
    object = Seurat_subset,
    group.by = "cluster_ann",
    macs2.path = "/users/csoto/.conda/envs/macs2_conda3_env/bin/macs2",
    verbose = TRUE
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


##==============================================================================

# Get current cluster identities
current_idents <- as.character(Idents(seurat_obj))

# Replace LHb cluster names with 'LHb_merged'
merged_idents <- ifelse(current_idents %in% LHb_clusters_to_merge, 
                        "LHb_merged", 
                        current_idents)

# Assign new identities to the Seurat object
Idents(seurat_obj) <- merged_idents

##==============================================================================


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


