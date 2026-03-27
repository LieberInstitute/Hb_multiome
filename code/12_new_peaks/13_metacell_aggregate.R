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
graining_level = 75

dir.create(dirname(seur_out_path), showWarnings = FALSE)
dir.create(dirname(plot_path), showWarnings = FALSE)

seur = qs_read(seur_in_path)

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

qs_save(seur_meta, seur_out_path)

#   Show metacells on UMAP dimensions
pdf(plot_path)
DimPlotSC(
    seur, seur_meta, reduction = "wnn.umap", sc.col = "mid_cluster",
    metacell.col = "mid_cluster"
)
dev.off()

FetchData(seur_meta, c("mid_cluster", "mid_cluster_purity")) |>
    as_tibble() |>
    pull(mid_cluster_purity) |>
    summary()

################################################################################
#   Rebuild missing parts of the metacell-level object
################################################################################

Fragments(seur_meta) = Fragments(seur)

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
