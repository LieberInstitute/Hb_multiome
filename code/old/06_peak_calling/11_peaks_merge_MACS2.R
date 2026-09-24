########################################################################
## Workflow for merging the peaks and creating a new assay with merged peaks
## 
## Authors. CSC
## Date. Sep 03, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
## I use BSgenome.Hsapiens.UCSC.hg38 for extracting DNA motifs, k-mers, sequence-based features and compute Tn5 bias correction
library("BSgenome.Hsapiens.UCSC.hg38")  # full reference genome sequence
library("GenomicRanges") 
library("purrr")
library("tidyverse")
library("tidyr")
library("stringr")
library("here")


## ATAC function's helper used globally
source(here("code", "06_peak_calling", "atac_custom_functions", "atac_normalization_helpers.R"))
# ls()

#===============================================================================
# resolution_level = "Broad"    # 8 cell-types
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

p_met = "spearman"
w_size = "5e5"

# ## read input arguments
# args = commandArgs(trailingOnly = TRUE)
# resolution_level <- args[2]

## Fixed resolution for this dataset:
resolution_level = "Mid" 

if (length(resolution_level)) {
    message(
        "Processing job for resolution_level:\n",
        resolution_level
    )
    f_sufix <- paste0(".", p_met, ".", w_size, ".cells_filtered_2perc")
} else {
    message("Input argument missed")
    stop()
}


# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
inputRDS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "22_add_mid_level_clustering" # Seurat multiome final version with annotations
)
input_macs_file <- here(
    "processed-data",
    "06_peak_calling",
    "01_call_peaks_MACS2",
    "macs_peaks_Mid_resolution.csv"
)
input_genes_filtered_file <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_MACS2",
    "rna_filtered_genes_2perc_cells.csv"
)
output_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "11_peaks_merge_MACS2"
)

if (!dir.exists(output_Dir)) {
    dir.create(output_Dir)
}

##==============================================================================
## Load Seurat / macs peaks / filtered genes. And make verification
 
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
SeuratOBJ
DefaultAssay(SeuratOBJ) <- "RNA"

# colnames(SeuratOBJ@meta.data)
msg <- switch(resolution_level,
    "Fine" = paste0("Fine-level Cell-Types:\n ", paste(sort(unique(SeuratOBJ$cluster_ann)), collapse = "\n")),
    "Broad" = paste0("Broad-level Cell-Types:\n", paste(sort(unique(SeuratOBJ$merged_cluster)), collapse = "\n")),
    "Mid" = paste0("Mid-level Cell-Types:\n", paste(sort(unique(SeuratOBJ$mid_cluster)), collapse = "\n"))
)

message(msg)
message("Processing ", length(Cells(SeuratOBJ)), " cells")

## set resolution_level
meta_col <- case_when(
    resolution_level=="Fine" ~ "cluster_ann",
    resolution_level=="Broad" ~ "merged_cluster",
    resolution_level=="Mid" ~ "mid_cluster"
)
 
## set desired idents as current level
SeuratOBJ <-  global_set_idents_from_meta(SeuratOBJ, meta_col = meta_col)
levels(SeuratOBJ)
clusters <- levels(SeuratOBJ)
message("Seurat loaded and ready!")
clusters

##==============================================================================
## load filtered genes expressed at least in 2% of cells & peaks from CallPeaks()

if (file.exists(input_genes_filtered_file)) {
    keep_genes <- read.csv(input_genes_filtered_file)$x
    head(keep_genes)
    message("Loaded ", length(keep_genes), " pre-filtered genes ...")
} else {
    stop(paste("File not found:", input_genes_filtered_file))
}


##==============================================================================
## load macs2 peaks

message("Loading ", basename(input_macs_file))

if (file.exists(input_macs_file)) {
    peaks_df <- read.csv(input_macs_file)
    message("File loaded!")
} else {
    stop(paste("File not found:", input_macs_file))
}
colnames(peaks_df)

## verify shorter and longer peaks
quantile(peaks_df$width, probs = c(0.01, 0.02, 0.05, 0.25, 0.5, 0.75, 0.95, 0.98, 0.99, 0.995, 0.997))
# 1%   2%   5%  25%  50%  75%  95%  98%  99% 
# 200  200  205  260  362  545 1131 1515 1766 
peaks_df[peaks_df$width < 200, ]
# <0 rows> (or 0-length row.names)
nrow(peaks_df[peaks_df$width > 2000, ])
# 1756
## Avoid very small or very large peaks (< 20bp or > 10kb). Improve quality
x <- nrow(peaks_df)
peaks_df <- peaks_df[peaks_df$width >= 20 & peaks_df$width <= 10000, ]
y <- nrow(peaks_df)

message("Peaks removed ", (x - y))
message("Peaks kept: ", (x - (x-y)), " (", round((y * 100) / x, digits = 2), "%)")

## CSV has columns like: seqnames, start, end, etc
peaks_gr <- makeGRangesFromDataFrame(
    peaks_df,
    keep.extra.columns = TRUE,       # keep additional metadata columns
    seqnames.field    = "seqnames",  # adjust if column name differs
    start.field       = "start",
    end.field         = "end"
)

## Ensure your fragment files and peaks use the same UCSC-style naming style
seqlevelsStyle(peaks_gr) <- "UCSC"

# inspect
class(peaks_gr)  # [1] "GenomicRanges"
length(peaks_gr) # 355127 -> mid resolution
head(peaks_gr)

message("MACS2 peaks prepared!")


##==============================================================================
## - merge peaks 
## - quantify peaks from CallPeaks output fragments

message("Stating Quantification ...")

## Set ATAC assay 
Assays(SeuratOBJ)
DefaultAssay(SeuratOBJ) <- "ATAC"

## Use *all* fragment files: meaning from all the samples 
frags_list <- Fragments(SeuratOBJ)
length(frags_list) # 10

message("Merging peaks ...")

# mapping the new, unified peaks back to the original individual peaks. This allows you to retain a link to the original peak 
union_peaks <- GenomicRanges::reduce(peaks_gr, min.gapwidth = 101, with.revmap = TRUE) 

# fast confirmation
length(union_peaks) # [1] 351037
head(union_peaks)
revmap <- mcols(union_peaks)$revmap  # IntegerList: indices of contributing peaks

## =============================================================================
## log summary and stat verificacion

message("Original number fo peaks:")
p_total <- length(peaks_gr)            # 355127 - original count
message("Unified number fo peaks:", p_total)
p_unified <- length(union_peaks)         # 351037 - unified count (should be <= original)
message("Peak difference:", p_unified)
p_total - p_unified # 4099
message("Merged peaks %:")
round((p_total - p_unified) *  100 / p_total, digits = 2)

## verification
any(width(union_peaks) <= 0)  # should be FALSE
any_overlaps <- any(countOverlaps(peaks_gr, peaks_gr) > 1)
any_overlaps # FALSE
is_disjoint <- isDisjoint(peaks_gr, ignore.strand = TRUE) 
is_disjoint # [1] TRUE / meaning none of the peak intervals overlapped with each other
n_dups <- sum(duplicated(peaks_gr)) # 0 duplicate intervals
n_dups # 0 

stopifnot(length(frags_list) >= 1)
all_cells_frag <- unique(do.call(c, lapply(frags_list, Cells)))
length(all_cells_frag) # 55516
barcodes_in_frags <- unique(unlist(all_cells_frag, Cells(SeuratOBJ)))
mean(colnames(SeuratOBJ) %in% barcodes_in_frags)  # should be ~1.0

## =============================================================================

message("Quantifying merged peaks ...")

mat_macs2_peaks <- FeatureMatrix(
    fragments = frags_list,              # list of Fragment objects
    features  = union_peaks,             # GRanges from CallPeaks
    cells     = colnames(SeuratOBJ),     # all barcodes across samples
    verbose   = TRUE
)

f_name <- paste0("mtx_merged_peaks_cell_level_", resolution_level, "_resolution.rds")
rds_name <- here(output_Dir, f_name)
saveRDS(mat_macs2_peaks, file = rds_name)
# To load it back:
# mat_macs2_peaks <- readRDS(rds_name)

message("ATAC-Counts saved ...")

# Verification
identical(colnames(mat_macs2_peaks), colnames(SeuratOBJ))
head(mat_macs2_peaks)



##==============================================================================
## Redo new ATAC assay with merged peaks

atac_macs2 <- CreateChromatinAssay(
    counts     = mat_macs2_peaks,
    ranges     = union_peaks,
    annotation = tryCatch(Annotation(SeuratOBJ), error = function(e) NULL)
)
atac_macs2

## add the new chromatin object to existing Seurat
SeuratOBJ[["ATAC_macs2_merged"]] <- atac_macs2
Assays(SeuratOBJ)

message("Seurat with new chromatin assay done: ATAC_macs2_merged")

DefaultAssay(SeuratOBJ) <- "ATAC_macs2_merged"

message("Starting GC content correction ... ")

# Set genome
genome <- BSgenome.Hsapiens.UCSC.hg38

SeuratOBJ <- RegionStats(
    object = SeuratOBJ,
    assay = "ATAC_macs2_merged",
    genome = genome
)

# Filter cells / QC ?
# before_peaks <- sum(SeuratOBJ[["nCount_ATAC"]])
# tmp_Sobj <- subset(SeuratOBJ, subset = nCount_ATAC < 75000 & nCount_ATAC > 1000)
# after_peaks <- sum(SeuratOBJ[["nCount_ATAC"]])

SeuratOBJ <- global_rebuild_atac_normalization(SeuratOBJ, "ATAC_macs2_merged")

message("GC content correction and normalization done!")

## Milestone: save whole Seurat for further analysis
f_name <- paste0("Seurat_peaks_merged_cell_level_", resolution_level, "_resolution.rds")
rds_name <- here(output_Dir, f_name)
saveRDS(SeuratOBJ, file = rds_name)

message("Seurat with new merged_peaks (chromatin assay) added and saved!")

message("All done!!!")


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

