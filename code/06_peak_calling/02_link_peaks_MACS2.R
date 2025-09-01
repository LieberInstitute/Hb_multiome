########################################################################
## Compute LinkPeaks() from CallPeaks() - MACS2
## Count macs2 peaks, make GC bias corrections and normalize new chromatin assay
## Next compute LinkPeaks correlations by cluster (multiome WNN)
## Save genes-filtered, peaks-mtx and update Seurat
## 
## Authors. CSC
## Date. August 18, 2025
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
# resolution_level = "Fine"     # 42 clusters (small clusters - not run)
# resolution_level = "Broad"    # 8 cell-types
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

p_met = "spearman"
w_size = "5e5"

## read input arguments
args = commandArgs(trailingOnly = TRUE)
resolution_level <- args[2]

## for testing:
# resolution_level = "Mid" 

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
    "01_call_peaks_MACS2"
)
output_cvsDir <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_MACS2" # local peaks redo with Signac::CallPeaks()
)
output_RDS <- here(
    "processed-data",
    "06_peak_calling",
    "02_link_peaks_MACS2"
)

if (!dir.exists(output_cvsDir)) {
    dir.create(output_cvsDir)
}
if (!dir.exists(output_RDS)) {
    dir.create(output_RDS)
}

##==============================================================================
## Load Seurat / macs peaks / filtered genes. And make verification

## macs2 peaks
macs2_peaks_csv <- here(input_macs_file, paste0("macs_peaks_Mid_resolution.csv"))

## Seurat with final ct - Visium HD corrected 
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
## filter genes to those expressed in 2% of cells

rna_counts <- GetAssayData(SeuratOBJ, assay="RNA", layer="data")
length(rownames(rna_counts)) # [1] 36601
#length(rownames(rna_counts)[Matrix::rowSums(rna_counts > 0)]) # 34738
#length(rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02]) # 34738

# filter rna count expressed in at least 2% of the cells
keep_genes <- rownames(rna_counts)[Matrix::rowSums(rna_counts > 0) > 0.02 * ncol(rna_counts)]
length(keep_genes) # in count: [1] 15896
# Note. By default behavior on Seurat v5, after remove cells, ATAC assay was removed, which will be attached later
#SeuratOBJ_subset <- SeuratOBJ
#SeuratOBJ_subset <- SeuratOBJ_subset[keep_genes, ]  
# Convert GRanges to data frame and save for further analysis
head(keep_genes)

## save genes to avoid redundant step on subsequent scripts
f_name <- paste0("rna_filtered_genes_2perc_cells.csv")
write.csv(
    keep_genes, 
    here(output_cvsDir, f_name),
    row.names = FALSE
)

message("Filtered gene names saved ...")


##==============================================================================

# load macs2 peaks

message("Loading ", basename(macs2_peaks_csv))

if (file.exists(macs2_peaks_csv)) {
    peaks_df <- read.csv(macs2_peaks_csv)
    message("File loaded!")
} else {
    stop(paste("File not found:", macs2_peaks_csv))
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
## Avoid very small or very large peaks (< 20bp or > 2kb). Improve quality
x <- nrow(peaks_df)
peaks_df <- peaks_df[peaks_df$width >= 200 & peaks_df$width <= 2000, ]
y <- nrow(peaks_df)

message("Peaks removed ", (x - y))
message("Peaks kept: ", (x - (x-y)), " (", round((y * 100) / x, digits = 2), "%)")

# CSV has columns like: seqnames, start, end, etc
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
head(peaks_gr)
length(peaks_gr) # 355127 -> mid resolution
head(peaks_gr)
# GRanges object with 6 ranges and 1 metadata column:
# seqnames        ranges strand |         peak_called_in
# <Rle>     <IRanges>  <Rle> |            <character>
# [1]     chr1 181329-181534      * |             Inhib.Thal

message("MACS2 peaks prepared!")


##==============================================================================

## Pre-processing: create new ATAC object from CallPeaks output fragments

message("Stating Quantification ...")

## Set ATAC assay 
Assays(SeuratOBJ)
DefaultAssay(SeuratOBJ) <- "ATAC"

# Use *all* fragment files: meaning from all the samples 
## pick-up smallest cell-type size to run test:
# tmp <- purrr::map(seurat_subsets, ~ length(Cells(.x)))
# tmp
# # get the name of the object with minimum number of cells
# names(tmp)[which.min(tmp)]
# small_idx_to_test <- which.min(tmp)
# as.integer(small_idx_to_test)

frags_list <- Fragments(SeuratOBJ)
length(frags_list) # 10

## verification
stopifnot(length(frags_list) >= 1)
all_cells_frag <- unique(do.call(c, lapply(frags_list, Cells)))
length(all_cells_frag) # 55516
barcodes_in_frags <- unique(unlist(all_cells_frag, Cells(SeuratOBJ)))
mean(colnames(SeuratOBJ) %in% barcodes_in_frags)  # should be ~1.0

## Quantify peaks detected by Signac::CallPeaks() - not unified peaks
mat_macs2_peaks <- FeatureMatrix(
    fragments = frags_list,              # list of Fragment objects
    features  = peaks_gr,                # GRanges from CallPeaks
    cells     = colnames(SeuratOBJ),     # all barcodes across samples
    verbose   = TRUE
)

# ## if crashes for memory, I will try to filter for peaks occurring in a minimum number of cells
# peak_counts <- Matrix::rowSums(mat_macs2_peaks > 0)
# keep_peaks <- peak_counts > 50 
# mat_macs2_peaks <- mat_macs2_peaks[keep_peaks, ]
# or
# Avoid very small or very large peaks (< 20bp or > 2kb)

message("Peaks quantification done!")

f_name <- paste0("mtx_peaks_cell_level_", resolution_level, "_resolution.rds")
rds_name <- here(output_RDS, f_name)
saveRDS(mat_macs2_peaks, file = rds_name)
# To load it back:
# mat_macs2_peaks <- readRDS(rds_name)

message("Chromatin counts saved ...")

# identical(colnames(mat_macs2_peaks), colnames(SeuratOBJ))
# head(mat_macs2_peaks)
# # > head(mat_macs2_peaks) # ge. Mid_level: 18 x 55,516 cells
# # 6 x 55516 sparse Matrix of class "dgCMatrix"
# # [[ suppressing 34 column names ‘S04_AAACAGCCAGAATGAC-1’, ‘S04_AAACAGCCAGCAAGGC-1’, ‘S04_AAACATGCACCTGGTG-1’ ... ]]

# create new ATAC assay with macs2 Granges
atac_macs2 <- CreateChromatinAssay(
    counts     = mat_macs2_peaks,
    ranges     = peaks_gr,
    annotation = tryCatch(Annotation(SeuratOBJ), error = function(e) NULL)
)
atac_macs2

## add the new chromatin object to existing Seurat 
SeuratOBJ[["ATAC_macs2"]] <- atac_macs2
Assays(SeuratOBJ)
# [1] "RNA"        "ATAC"       "ATAC_macs2"

message("Seurat with new chromatin assay: ATAC_macs2")

DefaultAssay(SeuratOBJ) <- "ATAC_macs2"

# Make GC content correction and run TF-IDF normalization, runSVD & 

message("Starting GC content correction ... ")

# Set genome
genome <- BSgenome.Hsapiens.UCSC.hg38

# Compute GC content for each peak
SeuratOBJ <- RegionStats(
    object = SeuratOBJ,
    assay = "ATAC_macs2",
    genome = genome
)

message("GC content correction and normalization done!")

# Filter cells / QC ? 

# before_peaks <- sum(SeuratOBJ[["nCount_ATAC"]])
# tmp_Sobj <- subset(SeuratOBJ, subset = nCount_ATAC < 75000 & nCount_ATAC > 1000)
# after_peaks <- sum(SeuratOBJ[["nCount_ATAC"]])

SeuratOBJ <- global_rebuild_atac_normalization(SeuratOBJ, "ATAC_macs2")

## milestone
f_name <- paste0("Seurat_peaks_macs2_cell_level_", resolution_level, "_resolution.rds")
rds_name <- here(output_RDS, f_name)
saveRDS(SeuratOBJ, file = rds_name)
#readRDS(SeuratOBJ, file = rds_name)

message("Seurat with new chromatin assay GC bias corrected and normalized saved!")


##==============================================================================
## Compute cell-type specific (local) Link peak-genes

## find peaks by cluster correlated with the expression of nearby genes
# access each subset by name

message("Making list of Seurat subsets:", resolution_level)

seurat_subsets <- list()

for (clust in clusters) {
    message("Subsetting cluster: ", clust)
    
    seurat_subsets[[clust]] <- subset(
        SeuratOBJ,
        idents = clust
    )
}

message("Seurat subsets by cell-type arranged: ", length(seurat_subsets))

message("Starting LinkPeaks by cluster ... ")

set.seed(22082025)

for (seurat_cluster in names(seurat_subsets)) {
    # seurat_cluster = "Endo"

    print(paste("Processing seurat cluster: ", seurat_cluster))

    seurat_subset <- seurat_subsets[[seurat_cluster]]
    seurat_subset

    print(paste("Processing Peaks for ", seurat_cluster, "\n",
                "Total cells found:", length(Cells(seurat_subset))))

    atac <- LinkPeaks(
        object = seurat_subset,
        peak.assay = "ATAC_macs2",
        expression.assay = "RNA",
        genes.use = keep_genes,
        method = p_met,
        distance = as.numeric(w_size)             # Only consider peaks within x kb of gene TSS
    )

    # This is for embedding the peak-gene links in the Seurat object directly, in case I need it
    # Links(SeuratOBJ[["ATAC_macs2"]]) <- Links(atac)
    
    print("Local link-peaks correlations completed!")

    ## inspect data
    print(head(Links(atac), n=3))

    ## prepare data to save cvs
    link_df <- as.data.frame(Links(atac))
    print(summary(link_df$score))

    f_name <- paste0(resolution_level, "_", seurat_cluster, "_local_link_peak_genes", f_sufix, ".csv")
    write.csv(
        link_df,
        file = here(output_cvsDir, f_name),
        row.names = FALSE
    )
    
    message("LinkPeaks saved: ", f_name)

}


message("MACS2 link-peaks correlations completed!")


message("All done!!!")


# library("slurmjobs")
# job_single(
#   "00_link_peaks",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
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
