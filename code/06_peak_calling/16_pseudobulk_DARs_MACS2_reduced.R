########################################################################
## Differential Accessibility (DA) using Seurat/Signac from pseudobulk assays on meged peaks
## - Peaks were computed with MACS2
##
## Authors. CSC
## Date. Sep 11, 2025
## Recommended resources mem=30GB
########################################################################

# library("EnsDb.Hsapiens.v86")  # Gene annotation (GTF-style), gene names, positions, TSSs, chr locations, etc.
# library("ggplot2")
# library("patchwork")
# library("tidyverse")
# library("stringr")
library("Seurat")
library("Signac")
library("dplyr")
# library("scales")
library("here")


# Testing spearman at 5e4 on macs2 peaks 
resolution_level = "Mid"
p_met = "spearman"
w_size = "5e5"
# FDR_thresh = 0.2
# score_thresh = 0.2

if (length(resolution_level)) {
    message("Processing job for peak-method:\n",
            p_met,
            "\nWindow-size\n",
            w_size)
} else {
    message("Input arguments missed")
    stop()
}

# Check/create directories
inputRDS_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "12_pseudobulk_MACS2"
)
output_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "16_pseudobulk_DARs_MACS2_reduced.R"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "16_pseudobulk_DARs_MACS2_reduced.R"
)


## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(output_Dir)) {
    dir.create(output_Dir)
}

# List all files matching the specific clustering resolution level
lst_peak_files <- list.files(
    path = inputRDS_Dir,
    pattern = "*.rds",
)
#lst_peak_files = list.files(path = input_cvsDir)
message("Link peak-genes files found:")
lst_peak_files
# [1] "mtx_merged_peaks_cell_level_Mid_resolution.rds"   
# [2] "Seurat_peaks_merged_cell_level_Mid_resolution.rds"


########################################################################
## Differential Accessibility Analysis in Seurat/Signac
## Input: Seurat object with pseudobulk RNA+ATAC assay and reduced peaks
########################################################################

##==============================================================================
## Load Seurat / macs peaks / filtered genes. And make verification

Seurat_base_name <- "Mid_pseudobulk.spearman.5e5_merged_peaks.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ_pb <- readRDS(here(inputRDS_Dir, Seurat_base_name))
# check
SeuratOBJ_pb
# An object of class Seurat 
# 366933 features across 169 samples within 2 assays 
# Active assay: ATAC_macs2_merged_pseudo (351037 features, 333643 variable features)
# 2 layers present: counts, data
# 1 other assay present: RNA
# 1 dimensional reduction calculated: lsi

colnames(SeuratOBJ_pb@meta.data)
# [1] "orig.ident"                        "nCount_RNA"                       
# [3] "nFeature_RNA"                      "nCount_ATAC_macs2_merged_pseudo"  
# [5] "nFeature_ATAC_macs2_merged_pseudo"
head(SeuratOBJ_pb@meta.data)


levels(SeuratOBJ_pb)
unique(Idents(SeuratOBJ_pb))

message("Pseudobulk groups:")
unique(SeuratOBJ_pb[["orig.ident"]])

## Define assay and metadata grouping
PSEUDO_ATAC_ASSAY <- "ATAC_macs2_merged_pseudo"  
grouping_var <- "orig.ident"    

# Set default assay
DefaultAssay(SeuratOBJ_pb) <- PSEUDO_ATAC_ASSAY


## Find differentially accessible peaks between cell-types
## Loop over multiple identities

# testing 
# da_res <- FindMarkers(
#     object = SeuratOBJ_pb,
#     ident.1 = "LHb.1",
#     ident.2 = NULL,              # NULL = rest of cells
#     only.pos = FALSE,            # get both up and down
#     min.pct = 0.1,               # require peak present in 10% of cells
#     test.use = "LR"              # likelihood ratio test (works well for binary peak presence)
#     # latent.vars = "nCount_ATAC"  # correct for sequencing depth if needed => this is pb data
# )


## provide as input reduced peaks, so, this set should be consistent across cells (important for DA)

cluster_ids <- levels(SeuratOBJ_pb)

message("Testing")
cluster_ids

all_da <- lapply(cluster_ids, function(ct) {
    res <- FindMarkers(
        seurat_obj,
        ident.1 = ct,
        only.pos = TRUE,
        test.use = "LR"              # recommended for binary peak data in Signac / Opt. poisson and negbinom
    )
    res$cluster <- ct
    res$peak <- rownames(res)
    res$FDR <- p.adjust(res$p_val, method = "BH")
    return(res)
})

da_results <- bind_rows(all_da)
write.csv(da_results, "DA_all_clusters.csv", row.names = FALSE)

# #------------------------------------------------------
# # 5. Optional: Annotate peaks
# #------------------------------------------------------
# # Example: distance to nearest gene TSS
# annotations <- ClosestFeature(seurat_obj, regions = rownames(seurat_obj))
# da_results_annot <- left_join(da_results, annotations, by = c("peak" = "query_region"))
# write.csv(da_results_annot, "DA_all_clusters_annotated.csv", row.names = FALSE)
