########################################################################
## Call LinkPeaks() on pseudobulk at Mid cell-type level
## Parse Seurat by ct due memory issues
## 
## Authors. CSC
## Date. Sep 07, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("BSgenome.Hsapiens.UCSC.hg38")
library("GenomicRanges") 
library("purrr")
library("tidyverse")
library("tidyr")
library("stringr")
library("here")


## ATAC function's helper used globally
source(here("code", "06_peak_calling", "multiome_custom_functions", "multiome_idents_normalization_helper.R"))
# ls()

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

p_met = "spearman"
w_size_num <- 5e5 
w_size = "5e5" # for filenames only 
resolution_level = "Mid" 
#cluster_name <- "Endo"

## read input arguments
args = commandArgs(trailingOnly = TRUE)
cluster_name <- args[2]
# cluster_name = "Inhib.Thal"
if (is.na(cluster_name) || !nzchar(cluster_name)) stop("Missing cluster_name argument")

## Check/create directories
inputRDS_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "12_pseudobulk_MACS2"
)
output_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "13_pseudobulk_LinkPeaks_MACS2_split_ct"
)

if (!dir.exists(output_Dir)) {
    dir.create(output_Dir)
}


##==============================================================================
## Load Seurat / macs peaks / filtered genes. And make verification

SeuratOBJ <- readRDS(here(inputRDS_Dir, "Mid_pseudobulk.spearman.5e5.rds"))
pb_obj <- SeuratOBJ
pb_obj
# An object of class Seurat 
# 371023 features across 169 samples within 2 assays 
# Active assay: ATAC_macs2_pseudo (355127 features, 337546 variable features)
# 2 layers present: counts, data
# 1 other assay present: RNA
# 1 dimensional reduction calculated: lsi

DefaultAssay(pb_obj) <- "RNA"
levels(pb_obj)
# [1] "Astrocyte"  "Endo"       "Excit.Thal" "Inhib.Thal" "LHb.1"     
# [6] "LHb.1.3"    "LHb.1.3.4"  "LHb.2.7"    "LHb.4"      "LHb.7"     
# [11] "MHb.1"      "MHb.1.2"    "MHb.2"      "MHb.3"      "Microglia" 
# [16] "Oligo"      "OPC"        "Thal" 

message("Assays present: ", paste(Assays(pb_obj), collapse=", "))
# Assays present: RNA, ATAC_macs2_pseudo


##==============================================================================
# genome for RegionStats
genome <- BSgenome.Hsapiens.UCSC.hg38

# derive cluster labels
if (!"cluster_pb" %in% colnames(pb_obj@meta.data)) {
    pb_obj$cluster_pb <- sub("[-_].*$", "", colnames(pb_obj))
}
pb_obj$cluster_pb <- trimws(pb_obj$cluster_pb)
Idents(pb_obj) <- pb_obj$cluster_pb
clusters_pb <- levels(Idents(pb_obj))

message(length(clusters_pb), " Clusters in pseudobulk: ", paste(clusters_pb, collapse=", "))

stopifnot(cluster_name %in% clusters_pb)

# test all genes in pb_obj
genes_for_lp <- rownames(pb_obj[["RNA"]])


##=====================================================================

# Use the pseudobulk ATAC assay that exists
atac_assay_name <- if ("ATAC_macs2_pseudo" %in% Assays(pb_obj)) "ATAC_macs2_pseudo" else stop("ATAC_macs2_pseudo does not exist!")
DefaultAssay(pb_obj) <- atac_assay_name

# peaks from assay and peak->cluster map (from metadata)
peaks_gr <- granges(pb_obj[[atac_assay_name]])
head(peaks_gr)
# GRanges object with 6 ranges and 1 metadata column:
#     seqnames        ranges strand |         peak_called_in
# <Rle>     <IRanges>  <Rle> |            <character>
# [1]     chr1 181329-181534      * |             Inhib.Thal
# [2]     chr1 191217-191619      * | OPC,Oligo,Inhib.Thal..
# [3]     chr1 629146-629354      * |              Astrocyte
# [4]     chr1 629811-630032      * | LHb.7,Astrocyte,Olig..
# [5]     chr1 630189-630389      * |             Oligo,Endo
# [6]     chr1 632189-632410      * |              Astrocyte
peaks_gr$peak_id <- Signac::GRangesToString(peaks_gr)

# MACS2 GRanges stored as metadata column 'peak_called_in' in peaks_gr
stopifnot("peak_called_in" %in% colnames(mcols(peaks_gr)))

peaks_map <- as.data.frame(peaks_gr) |>
    transmute(peak_id = peak_id,
              peak_called_in = as.character(peak_called_in)) |>
    mutate(peak_called_in = str_split(peak_called_in, "\\s*,\\s*")) |>
    tidyr::unnest(peak_called_in) |>
    mutate(peak_called_in = trimws(peak_called_in)) |>
    filter(!is.na(peak_called_in), peak_called_in != "") |>
    distinct()

head(peaks_map)
# # A tibble: 6 × 2
# peak_id            peak_called_in
# <chr>              <chr>         
# 1 chr1-181329-181534 Inhib.Thal    
# 2 chr1-191217-191619 OPC           
# 3 chr1-191217-191619 Oligo         
# 4 chr1-191217-191619 Inhib.Thal    
# 5 chr1-191217-191619 LHb.4         
# 6 chr1-191217-191619 MHb.2 

peaks_by_cluster <- function(cluster) {
    peaks_map |> filter(peak_called_in == cluster) |> pull(peak_id) |> unique()
}

# subset pseudobulk to the requested cluster
seurat_subset <- subset(pb_obj, idents = cluster_name)
n_samples <- ncol(seurat_subset)
message("Processing ", cluster_name, " (", n_samples, " pseudobulk samples)")

if (n_samples < 3) stop("Too few samples for cluster ", cluster_name, "")


##==============================================================================

# keep only this cluster’s MACS2-called peaks present in assay
cluster_peaks <- intersect(peaks_by_cluster(cluster_name),
                           rownames(seurat_subset[[atac_assay_name]]))
if (length(cluster_peaks) == 0) stop("No MACS2 peaks for ", cluster_name, " in assay.")

# drop ultra-sparse peaks in this cluster (improves stability)
counts_mat <- GetAssayData(
    seurat_subset, 
    assay = atac_assay_name, 
    layer = "counts")[cluster_peaks, , drop = FALSE]

min_cells_sub <- max(3, floor(0.05 * n_samples))   # ≥5% or at least 3
keep_peaks_sub <- rownames(counts_mat)[Matrix::rowSums(counts_mat > 0) >= min_cells_sub]
if (length(keep_peaks_sub) == 0) stop("No peaks pass support filter in ", cluster_name, ".")


# Rebuild assay with kept peaks (preserves annotation cleanly)
peak_ranges <- Signac::StringToGRanges(keep_peaks_sub)

new_assay  <- CreateChromatinAssay(
    counts = counts_mat[keep_peaks_sub, , drop = FALSE],
    ranges = peak_ranges,
    annotation = tryCatch(Annotation(pb_obj[[atac_assay_name]]), error = function(e) NULL)
)
# ChromatinAssay data with 5225 features for 10 cells
# Variable features: 0 
# Genome: 
#     Annotation present: TRUE 
# Motifs present: FALSE 
# Fragment files: 0 

seurat_subset[[atac_assay_name]] <- new_assay
# Typically when assigning a new ChromatinAssay, a warning tells the new assay doesn’t have exactly the same feature set as the old one (ATAC_macs2_pseudo). That’s expected after merged or filtered peaks. Harmless.

DefaultAssay(seurat_subset) <- atac_assay_name

# RegionStats (attach GC%, width) — needed for bias correction
seurat_subset <- RegionStats(
    object = seurat_subset,
    assay  = atac_assay_name,
    genome = genome
)
# occasionally the peaks set contains seqnames not present in the genome package (e.g., random contigs, chrM if filtered, unplaced scaffolds). Those peaks simply won’t get stats assigned. Harmless

seurat_subset <- global_rebuild_atac_normalization(seurat_subset, atac_assay_name)
# Seurat auto-detects that only 9 SVD are valid (effective rank of the pb count matrix)
# silently trims to 9 as I have only ~10 pseudobulk groups, and I can’t extract 50 components anyway.

# Genes present in this subset’s RNA
genes_sub <- intersect(genes_for_lp, rownames(seurat_subset[["RNA"]]))
if (length(genes_sub) == 0) stop("No RNA genes to test in subset.")

message("Summary:",
        "\n  Cell-type = ", cluster_name,
        "\n  Peaks = ", length(keep_peaks_sub),
        "\n  Genes = ", length(genes_sub),
        "\n  min.cells = ", min_cells_sub,
        "\n  distance = ", w_size_num)
seurat_subset
# Summary:
#     Cell-type = LHb.1.3
# Peaks = 3409
# Genes = 15896
# min.cells = 3
# distance = 5e+05
# An object of class Seurat 

## ============================================================================/

message("Computing peak-gene correlations on pseudobulk by cell-type ")

DefaultAssay(seurat_subset) <- "RNA"  # not strictly required but conventional

# retain all links for downstream multiple testing correction (FDR/HB)
seurat_subset <- LinkPeaks(
    object           = seurat_subset,
    peak.assay       = atac_assay_name,
    expression.assay = "RNA",
    genes.use        = genes_sub,
    distance         = w_size_num,
    min.cells        = min_cells_sub,
    pvalue_cutoff    = 1,   # keep everything / keep only links with pvalue <= pvalue.cutoff
    score_cutoff     = 0,   # keep both positive and negative scores
    method           = p_met
)

# Extract links and compute FDR
lk_gr <- Links(seurat_subset[[atac_assay_name]])
lk_df <- as.data.frame(lk_gr)

# Add FDR column
lk_df$FDR <- p.adjust(lk_df$pvalue, method = "BH")
lk_df$cluster <- cluster_name

message("Links found: ", nrow(lk_df))
# Links found: 28251
if (nrow(lk_df) > 0) {
    print(head(lk_df[, c("peak", "gene", "score", "pvalue", "FDR")], 5))
}

# Save as CSV
out_csv <- here(output_Dir, paste0(resolution_level, "_", cluster_name,
                                   "_pseudobulk_link_peak_genes.csv"))
write.csv(lk_df, out_csv, row.names = FALSE)
message("Links: ", nrow(lk_df), " (saved: ", basename(out_csv), ")")

# Save subset if desired
f_name <- paste0(resolution_level, "_", cluster_name,
                 "_pseudobulk_seurat_subset.rds")
saveRDS(seurat_subset, file = here(output_Dir, f_name))

message("All done!!!")


## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

