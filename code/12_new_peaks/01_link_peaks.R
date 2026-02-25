library(getopt)
library(sessioninfo)
library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges) 
library(tidyverse)
library(here)
library(Matrix)
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
    "processed-data", "06_peak_calling", "12_pseudobulk_MACS2",
    "Mid_pseudobulk.spearman.5e5.rds"
)
out_path = here(
    "processed-data", "12_new_peaks", "01_link_peaks",
    sprintf("%s_%s.csv.gz", opt$target_cell_type, opt$other_cell_type)
)
atac_assay = "ATAC_macs2_pseudo"
rna_assay = "RNA"

message(Sys.time(), ' | Loading Seurat object')
seur = readRDS(seur_path)

message(Sys.time(), ' | Rebuilding ATAC assay with subset peaks')

#   Set Idents and subset to the other cell type
seur$cluster_pb = sub("[-_].*$", "", colnames(seur))
seur$cluster_pb = trimws(seur$cluster_pb)
Idents(seur) = seur$cluster_pb
seur = subset(seur, idents = opt$other_cell_type)

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

#   Filter peaks to those present in at least 5% of cells in the other cell type
counts_mat = GetAssayData(seur, assay = atac_assay, layer = "counts")[
    peaks_called, , drop = FALSE
]
min_cells = max(3, floor(0.05 * ncol(seur)))   # ≥5% or at least 3
keep_peaks = rownames(counts_mat)[
    (Matrix::rowSums(counts_mat > 0) >= min_cells) &
    (sparseMatrixStats::rowSds(counts_mat) > 0)
]
stopifnot(length(keep_peaks) > 0)

# Rebuild assay with kept peaks (preserves annotation cleanly)
new_assay = CreateChromatinAssay(
    counts = counts_mat[keep_peaks, , drop = FALSE],
    ranges = StringToGRanges(keep_peaks),
    annotation = Annotation(seur[[atac_assay]])
)
seur[[atac_assay]] = new_assay
DefaultAssay(seur) = atac_assay

# RegionStats (attach GC%, width) — needed for bias correction
message(Sys.time(), ' | Running RegionStats')
seur = RegionStats(
    object = seur, assay = atac_assay, genome = BSgenome.Hsapiens.UCSC.hg38
)

#   Re-normalize following Cynthia's code
message(Sys.time(), ' | Re-normalizing ATAC assay')
seur = RunTFIDF(seur, assay = atac_assay, method = 1, scale.factor = 10000)
seur = FindTopFeatures(
    seur, assay = atac_assay, min.cutoff = 'q5', verbose = TRUE
)
seur = RunSVD(seur, assay = atac_assay)

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
