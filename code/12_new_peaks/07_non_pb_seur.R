#   The new linked peaks in this directory were created from pseudobulked,
#   non-merged peaks. However for TRIPOD, we need non-pseudobulked, non-merged
#   peaks. A Seurat object of that sort was created at one point based on
#   investigating earlier scripts, but is currently missing from the expected
#   path. I checked that the peaks themselves (which do exist as a data.frame)
#   match IDs with the new linked peaks, so in this script I'll take those peaks
#   and form a non-pseudobulked, non-merged Seurat object to be used for TRIPOD

library(tidyverse)
library(Seurat)
library(Signac)
library(here)
library(qs2)
library(sessioninfo)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(harmony)

peak_path = here(
    "processed-data", "06_peak_calling", "01_call_peaks_MACS2",
    "macs_peaks_Mid_resolution.csv"
)
seur_in_path = here(
    "processed-data", "05_Clustering_ARCr", "24_smaller_seurat",
    "general_purpose_seur.qs2"
)
seur_out_path = here(
    "processed-data", "12_new_peaks", "07_non_pb_seur",
    "non_pb_seur.qs2"
)
example_peaks_path = here(
    "processed-data", "12_new_peaks", "01_link_peaks",
    "Astrocyte_Astrocyte.csv.gz"
)
plot_path = here(
    "plots", "12_new_peaks", "07_non_pb_seur",
    "harmony_convergence.pdf"
)

dir.create(dirname(seur_out_path), showWarnings = FALSE)
dir.create(dirname(plot_path), showWarnings = FALSE)

################################################################################
#   Read in cell-level Seurat and called peaks
################################################################################

message(Sys.time(), " | Reading in Seurat object and called peaks...")
seur = qs_read(seur_in_path)
DefaultAssay(seur) = "ATAC"

#   Read in called peaks
peaks_gr = read_csv(peak_path, show_col_types = FALSE) |>
    makeGRangesFromDataFrame(
        keep.extra.columns = TRUE,
        seqnames.field = "seqnames",
        start.field = "start",
        end.field = "end"
    )
seqlevelsStyle(peaks_gr) = "UCSC"

################################################################################
#   Rebuild ATAC assay with those peaks
################################################################################

message(Sys.time(), " | Rebuilding ATAC assay with new peaks...")
peaks_mat = FeatureMatrix(
    fragments = Fragments(seur),
    features  = peaks_gr,
    cells = colnames(seur),
    verbose = TRUE
)
new_assay = CreateChromatinAssay(
    counts = peaks_mat,
    ranges = peaks_gr,
    fragments = Fragments(seur),
    annotation = Annotation(seur)
)
seur[['ATAC']] = new_assay
DefaultAssay(seur) = "ATAC"

#   As a sanity check, validate that some example peaks from the new linked
#   peaks all exist in this new assay
example_peaks_df = read_csv(example_peaks_path, show_col_types = FALSE)
stopifnot(all(example_peaks_df$peak %in% rownames(seur[['ATAC']])))

################################################################################
#   Normalize and add stats (follows Cynthia's workflow)
################################################################################

message(Sys.time(), " | Normalizing and adding stats...")
seur = RegionStats(
    object = seur, assay = "ATAC", genome = BSgenome.Hsapiens.UCSC.hg38
)

#   The choice to rerun TFIDF and SVD is a subtle one-- we're essentially
#   reframing the genomic ranges that define each peak (from those chosen by
#   CellRanger originally to those called by MACS2). In theory since reductions
#   are at the cell level and the underlying fragment counts are the same, the
#   LSI embedding should be fairly similar before and after. I'm opting to rerun
#   this step though, as well as Harmony, as I feel it's slightly more precise
seur = RunTFIDF(
    seur, assay = "ATAC", method = 1, scale.factor = 10000
)   
seur = FindTopFeatures(
    seur, assay = "ATAC", min.cutoff = 'q5', verbose = TRUE
)
seur = RunSVD(seur, assay = "ATAC")

#   Following Cynthia's code here (https://github.com/LieberInstitute/Hb_multiome/blob/7f749651cb989c7f072c239a704b7f22933d4b4b/code/03_pseudobulking/08_harmony_CR_ARCr.R#L211-L219)
seur = RunHarmony(
    seur, group.by.vars = "orig.ident",
    reduction.save = "integrated.lsi.harmony", assay.use = "ATAC",
    reduction.use= 'lsi', plot_convergence = TRUE, early_stop = TRUE,
    project.dim = FALSE
)

################################################################################
#   Save the Seurat object
################################################################################

message(Sys.time(), " | Saving Seurat object...")
qs_save(seur, seur_out_path)

session_info()
