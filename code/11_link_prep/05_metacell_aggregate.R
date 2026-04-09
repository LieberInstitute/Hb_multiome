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
library(BSgenome.Hsapiens.UCSC.hg38)

seur_in_path = here(
    'processed-data', '11_link_prep', '02_rebuild_atac_assay',
    'cell_level_seur.qs2'
)
seur_out_path = here(
    "processed-data", "11_link_prep", "05_metacell_aggregate",
    "meta_seur.qs2"
)
plot_path = here(
    "plots", "11_link_prep", "05_metacell_aggregate",
    "metacell_UMAP.pdf"
)
graining_level = 20

dir.create(dirname(seur_out_path), showWarnings = FALSE)
dir.create(dirname(plot_path), showWarnings = FALSE, recursive = TRUE)

seur = qs_read(seur_in_path)

################################################################################
#   Form metacells
################################################################################

#   Form metacells using just the RNA, since even the harmonized ATAC data has
#   a notable donor structure, and the goal is to avoid donor-driven linked
#   peaks downstream
seur_meta = SCimplify_for_Seurat(
    seur, assay = "RNA", reduction = list("integrated.harmony"), 
    dims = list(1:30), gamma = graining_level,
    label = "refined_mid_cluster"
)
stopifnot(all(seur_meta@meta.data$refined_mid_cluster_purity == 1))

message("Donor purity:")
summary(seur_meta@meta.data$orig.ident_purity)

message("Number of metacells per cell type:")
table(seur_meta@meta.data$refined_mid_cluster)

#   Show metacells on UMAP dimensions
pdf(plot_path)
DimPlotSC(
    seur, seur_meta, reduction = "wnn.umap", sc.col = "refined_mid_cluster",
    metacell.col = "refined_mid_cluster"
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

seur_meta = RegionStats(
    object = seur_meta, assay = "ATAC", genome = BSgenome.Hsapiens.UCSC.hg38
)

#   ATAC normalization. Note we're using method 2 for TFIDF, whereas for the
#   cell-level object we used method 1. The SuperCell manuscript specifically
#   suggests method 2 for metacells
DefaultAssay(seur_meta) = "ATAC"
seur_meta = seur_meta |>
    RunTFIDF(assay = "ATAC", method = 2) |>   
    FindTopFeatures(assay = "ATAC", min.cutoff = 'q5', verbose = TRUE) |>
    RunSVD(assay = "ATAC")

#   RNA normalization
DefaultAssay(seur_meta) = "RNA"
seur_meta = seur_meta |>
    NormalizeData() |>
    FindVariableFeatures(selection.method = "vst")

qs_save(seur_meta, seur_out_path)

session_info()
