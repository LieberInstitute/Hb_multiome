library(sessioninfo)
library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges) 
library(tidyverse)
library(here)
library(Matrix)
library(sparseMatrixStats)

seur_path = here(
    "processed-data", "06_peak_calling", "12_pseudobulk_MACS2",
    "Mid_pseudobulk.spearman.5e5.rds"
)
out_path = here(
    "processed-data", "12_new_peaks", "01_link_peaks", "global.csv.gz"
)
atac_assay = "ATAC_macs2_pseudo"
rna_assay = "RNA"

message(Sys.time(), ' | Loading Seurat object')
seur = readRDS(seur_path)

message(Sys.time(), ' | Rebuilding ATAC assay with subset peaks')

#   Set Idents
seur$cluster_pb = sub("[-_].*$", "", colnames(seur))
seur$cluster_pb = trimws(seur$cluster_pb)
Idents(seur) = seur$cluster_pb

#   Filter peaks to those present in at least 5% of cells
counts_mat = GetAssayData(seur, assay = atac_assay, layer = "counts")
min_cells = floor(0.05 * ncol(seur)) # ≥5% 
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
    mutate(FDR = p.adjust(pvalue, method = "BH")) |>
    write_csv(out_path)

message("Memory usage:")
gc()

session_info()

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/
