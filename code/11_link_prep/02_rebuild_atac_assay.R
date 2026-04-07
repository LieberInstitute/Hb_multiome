#    After calling peaks, we need to reconstruct the ATAC assay to measure
#    accessibility in these peaks. Rebuild the cell-level Seurat object

library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(qs2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(harmony)

peak_path = here(
    'processed-data', '11_link_prep', '01_call_peaks',
    'macs3_peaks.csv.gz'
)
seur_in_path = here(
    'processed-data', '05_03_annotation_adjustments', '06_refined_annotations',
    'refined_annotation_multiomeHab_Seurat.qs2'
)
seur_out_path = here(
    'processed-data', '11_link_prep', '02_rebuild_atac_assay',
    'cell_level_seur.qs2'
)
plot_path = here(
    "plots", "11_link_prep", "02_rebuild_atac_assay",
    "harmony_convergence.pdf"
)

dir.create(dirname(plot_path), showWarnings = FALSE)
dir.create(dirname(seur_out_path), showWarnings = FALSE)

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

#   Drop 'scale.data' from the RNA assay; IMO we don't use it enough to justify
#   the amount of extra time it takes to load the object
seur[["RNA"]] = subset(
    seur[["RNA"]], cells = colnames(seur), layer = c("data", "counts")
)

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

################################################################################
#   Normalize and add stats (follows Cynthia's workflow)
################################################################################

message(Sys.time(), " | Normalizing and adding stats...")
seur = RegionStats(
    object = seur, assay = "ATAC", genome = BSgenome.Hsapiens.UCSC.hg38
)

#   The choice to rerun TFIDF and SVD is a subtle one-- we're essentially
#   reframing the genomic ranges that define each peak (from those chosen by
#   CellRanger originally to those called by MACS3). In theory since reductions
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
pdf(plot_path)
seur = RunHarmony(
    seur, group.by.vars = "orig.ident",
    reduction.save = "integrated.lsi.harmony", assay.use = "ATAC",
    reduction.use = 'lsi', plot_convergence = TRUE, project.dim = FALSE,
    early_stop = TRUE
)
dev.off()

################################################################################
#   Save the Seurat object
################################################################################

message(Sys.time(), " | Saving Seurat object...")
qs_save(seur, seur_out_path)

session_info()