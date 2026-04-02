library(getopt)
library(sessioninfo)
library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges) 
library(tidyverse)
library(here)
library(Matrix)
library(qs2)
library(sparseMatrixStats)

# Import command-line parameters
spec <- matrix(
    c(
        c("target_cell_type", "other_cell_type"),
        c("t", "o"),
        rep("1", 2),
        rep("character", 2),
        rep("Add variable description here", 2)
    ),
    ncol = 5
)
opt <- getopt(spec)

message("Using the following parameters:")
print(opt)

seur_path = here(
    "processed-data", "12_new_peaks", "13_metacell_aggregate",
    "seur_meta.qs2"
)
out_path = here(
    "processed-data", "12_new_peaks", "14_metacell_link_peaks",
    sprintf("%s_%s.csv.gz", opt$target_cell_type, opt$other_cell_type)
)
atac_assay = "ATAC"
rna_assay = "RNA"

dir.create(dirname(out_path), showWarnings = FALSE)

message(Sys.time(), ' | Loading Seurat object')
seur = qs_read(seur_path)

seur = RegionStats(
    object = seur, assay = "ATAC", genome = BSgenome.Hsapiens.UCSC.hg38
)

message(Sys.time(), ' | Subsetting ATAC assay with expressed peaks')

#   Subset to the other cell type
seur = subset(seur, subset = mid_cluster == opt$other_cell_type)

peaks_gr = granges(seur[[atac_assay]])
peaks_gr$peak_id = GRangesToString(peaks_gr)

#   Get peaks called in the target cell type
peaks_called = as.data.frame(peaks_gr) |>
    transmute(
        peak_id = peak_id,
        peak_called_in = as.character(peak_called_in)
    ) |>
    mutate(peak_called_in = str_split(peak_called_in, "\\s*,\\s*")) |>
    unnest(peak_called_in) |>
    mutate(peak_called_in = trimws(peak_called_in)) |>
    filter(!is.na(peak_called_in), peak_called_in != "") |>
    distinct() |>
    filter(peak_called_in == opt$target_cell_type) |>
    pull(peak_id)
stopifnot(all(peaks_called %in% rownames(seur[[atac_assay]])))

#   Filter peaks to those present in at least 5% of metacells in the other cell
#   type
counts_mat = GetAssayData(seur, assay = atac_assay, layer = "counts")[
    peaks_called, , drop = FALSE
]
min_cells = as.integer(0.05 * ncol(seur))
keep_peaks = rownames(counts_mat)[
    (Matrix::rowSums(counts_mat > 0) >= min_cells) &
    (sparseMatrixStats::rowSds(counts_mat) > 0)
]
stopifnot(length(keep_peaks) > 0)
message(
    sprintf(
        "Retained %.1f%% of peaks (%d / %d) based on expression",
        100 * length(keep_peaks) / length(peaks_called),
        length(keep_peaks),
        length(peaks_called)
    )
)
seur[['ATAC']] = subset(seur[[atac_assay]], features = keep_peaks)

#   Compute links without any filtering of outputs
message(Sys.time(), ' | Running LinkPeaks')
DefaultAssay(seur) = rna_assay
seur = LinkPeaks(
    object = seur, peak.assay = atac_assay, expression.assay = rna_assay,
    min.cells = min_cells, pvalue_cutoff = 1, score_cutoff = 0,
    method = "spearman"
)

#   Export linked peaks
message(Sys.time(), ' | Exporting to CSV')
Links(seur[[atac_assay]]) |>
    as.data.frame() |>
    as_tibble() |>
    mutate(
        FDR = p.adjust(pvalue, method = "BH"),
        target_cell_type = opt$target_cell_type,
        other_cell_type = opt$other_cell_type
    ) |>
    write_csv(out_path)

message("Memory usage:")
gc()

session_info()

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/
