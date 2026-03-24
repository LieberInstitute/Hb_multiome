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
atac_assay = 'ATAC'
graining_level = 75

dir.create(dirname(seur_out_path), showWarnings = FALSE)
dir.create(dirname(plot_path), showWarnings = FALSE)

seur = qs_read(seur_in_path)

#   Form metacells
seur_meta = SCimplify_for_Seurat(
    seur, assay = c("RNA", atac_assay),
    #   These are Harmony-corrected RNA PCs and ATAC LSI embeddings respectively
    reduction = list("integrated.harmony", "integrated.lsi.harmony"), 
    dims = list(rna.comp, adt.comp), gamma = graining_level,
    label = "mid_cluster"
)

#   Show metacells on UMAP dimensions
pdf(plot_path)
DimPlotSC(
    seur, seur_meta, reduction = "wnn.umap", sc.col = "mid_cluster",
    metacell.col = "mid_cluster"
)
dev.off()

qs_save(seur_meta, seur_out_path)

session_info()
