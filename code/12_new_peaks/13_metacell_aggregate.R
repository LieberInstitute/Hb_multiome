#   Essentially pseudobulk data by forming metacells using SuperCell 2.0.
#   Export the resulting Seurat object. Mostly follows the tutorial:
#   https://htmlpreview.github.io/?https://github.com/GfellerLab/SuperCell/blob/supercell-2.0/docs/tutorials/SuperCell2.0_PBMC_10x_multiome.html
#   except use the semi-supervised approach shown here:
#   https://htmlpreview.github.io/?https://github.com/GfellerLab/SuperCell/blob/supercell-2.0/docs/tutorials/SuperCell2.0_BM_CITE_seq.html#33_Semi-supervised_metacell_identifications

library(tidyverse)
library(Seurat)
library(Signac)
library(SuperCell)
library(here)
library(qs2)
library(sessioninfo)
library(GenomicRanges)

seur_in_path = here(
    "processed-data", "12_new_peaks", "07_non_pb_seur",
    "non_pb_seur.qs2"
)
seur_out_path = here(
    "processed-data", "12_new_peaks", "13_metacell_aggregate",
    "seur_meta.qs2"
)
plot_path = here(
    "plots", "12_new_peaks", "13_metacell_aggregate",
    "metacell_UMAP.pdf"
)
graining_level = 20

dir.create(dirname(seur_out_path), showWarnings = FALSE)
dir.create(dirname(plot_path), showWarnings = FALSE)

seur = qs_read(seur_in_path)

#   Manual cell-type edits while we decide the final definitions
cells_keep <- Cells(seur)[!seur@meta.data$mid_cluster %in% c("Thal", "LHb.7")]
seur <- subset(seur, cells = cells_keep)
seur@meta.data$mid_cluster = ifelse(
    seur@meta.data$mid_cluster %in% c("LHb.1", "LHb.1.3", "LHb.1.3.4"),
    "LHb.1.3.4",
    seur@meta.data$mid_cluster
)

################################################################################
#   Form metacells
################################################################################

#   Form metacells
seur_meta = SCimplify_for_Seurat(
    seur, assay = c("RNA", "ATAC"),
    #   These are Harmony-corrected RNA PCs and ATAC LSI embeddings respectively
    reduction = list("integrated.harmony", "integrated.lsi.harmony"), 
    dims = list(1:30, 2:30), gamma = graining_level,
    label = "mid_cluster"
)
stopifnot(all(seur_meta@meta.data$mid_cluster_purity == 1))

message("Donor purity:")
summary(seur_meta@meta.data$orig.ident_purity)

message("Number of metacells per cell type:")
table(seur_meta@meta.data$mid_cluster)

#   Show metacells on UMAP dimensions
pdf(plot_path)
DimPlotSC(
    seur, seur_meta, reduction = "wnn.umap", sc.col = "mid_cluster",
    metacell.col = "mid_cluster"
)
dev.off()

################################################################################
#   Rebuild missing parts of the metacell-level object
################################################################################

#   For now, we won't build aggregated fragments. See 
#   https://github.com/GfellerLab/SuperCell/issues/36

#   Retain peak metadata. This is for some reason quite complicated
gr_src <- granges(seur[["ATAC"]])
gr_tgt <- granges(seur_meta[["ATAC"]])
stopifnot(identical(rownames(seur[["ATAC"]]), rownames(seur_meta[["ATAC"]])))
stopifnot(identical(as.character(gr_src), as.character(gr_tgt)))
mcols(gr_tgt) <- mcols(gr_src)
methods::slot(seur_meta[["ATAC"]], "ranges") <- gr_tgt

#   ATAC normalization
DefaultAssay(seur_meta) = "ATAC"
seur_meta = seur_meta |>
    RunTFIDF(assay = "ATAC", method = 1, scale.factor = 10000) |>   
    FindTopFeatures(assay = "ATAC", min.cutoff = 'q5', verbose = TRUE) |>
    RunSVD(assay = "ATAC")

#   RNA normalization
DefaultAssay(seur_meta) = "RNA"
seur_meta = seur_meta |>
    NormalizeData() |>
    FindVariableFeatures(selection.method = "vst")

qs_save(seur_meta, seur_out_path)

session_info()
