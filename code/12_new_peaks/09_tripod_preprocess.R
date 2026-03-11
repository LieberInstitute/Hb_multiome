#   Preprocess data for TRIPOD, following many of the steps seen in the
#   vignette:
#   https://htmlpreview.github.io/?https://github.com/yuchaojiang/TRIPOD/blob/main/vignettes/preprocessing_pbmc.html

library(tidyverse)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(qs2)
library(here)
library(TRIPOD)
library(sessioninfo)

seur_path = here(
    "processed-data", "12_new_peaks", "08_chromVAR", "seur.qs2"
)
out_path = here(
    "processed-data", "12_new_peaks", "09_tripod_preprocess",
    "preprocessed_objects.qs2"
)

set.seed(0)
dir.create(dirname(out_path), showWarnings = FALSE)

seur = qs_read(seur_path)

tripod_seur = getObjectsForModelFit(object = seur, chr = paste0("chr", 1:22))

#   Here I manually reimplement TRIPOD::filterSeuratObject(), which hardcodes
#   the presence of an "SCT" assay. That isn't actually necessary for the TRIPOD
#   workflow, and we choose to normalize the RNA assay differently
genes = tripod_seur$transcripts.gr$gene_name
motifxTF = tripod_seur$motifxTF
seur@assays$RNA = subset(
    seur@assays$RNA, features = match(genes, rownames(seur@assays$RNA))
)
seur@assays$chromvar = subset(
    seur@assays$chromvar,
    features = match(motifxTF[, 1], rownames(seur@assays$chromvar))
)

#   Might want to manually reimplement this also, since it hardcodes many
#   specific choices about normalization and clustering (that are slightly
#   different than ours)
seur = processSeuratObject(
    object = seur, dim.rna = 1:50, dim.atac = 2:50, verbose = FALSE
)

################################################################################
#   Form metacells
################################################################################

cluster_df = optimizeResolution(
        object = seur, graph.name = "wsnn", assay.name = "WNN",
        resolutions = c(0.5, 1, 2, 4, 8, 16, 32), min.num = 20
    ) |>
    as_tibble()

#   The vignette doesn't provide an algorithmic way to select the best
#   resolution, but provides this recommendation:
#       "Our empirical approach is to select a resolution that gives 80 or more
#        metacells, the majority of which contain 20 or more single cells"
#   I essentially implement that logic algorithmically here
best_res = cluster_df |>
    dplyr::filter(num_clusters >= 80, num_below / num_clusters < 0.02) |>
    arrange(num_clusters) |>
    slice_head(n = 1) |>
    pull(resolution)

message(
    sprintf(
        "Selected resolution = %d. Here's the whole metric table:", best_res
    )
)
print(cluster_df)

seur = getClusters(
    object = seur, graph.name = "wsnn", algorithm = 3, resolution = best_res,
    verbose = FALSE
)

#   TRIPOD hardcodes reliance on pre-v5 Seurat objects. Do this for
#   compatibility
seur[['RNA']] = as(seur[['RNA']], "Assay")

metacell_seur = getMetacellMatrices(
    object = seur, cluster.name = "seurat_clusters"
)

################################################################################
#   Save in one list
################################################################################

pre_list = list(
    tripod_seur = tripod_seur,
    seur = seur,
    metacell_seur = metacell_seur
)

qs_save(pre_list, out_path)

session_info()
