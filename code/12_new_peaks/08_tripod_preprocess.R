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
    "processed-data", "12_new_peaks", "07_non_pb_seur",
    "non_pb_seur.qs2"
)
out_path = here(
    "processed-data", "12_new_peaks", "08_tripod_preprocess",
    "preprocessed_objects.qs2"
)

set.seed(0)

################################################################################
#   First filter to standard chromosomes (in ATAC)
################################################################################

seur = qs_read(seur_path)

DefaultAssay(seur) = "ATAC"
gr_std = keepStandardChromosomes(
    granges(seur[['ATAC']]), pruning.mode = "coarse"
)
peaks_keep = paste(seqnames(gr_std), start(gr_std), end(gr_std), sep = '-')
seur[['ATAC']] = subset(seur[['ATAC']], features = peaks_keep)

################################################################################
#   Run ChromVAR
################################################################################

pwm_set = getMatrixSet(
    x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE)
)
motif_mat = CreateMotifMatrix(
    features = granges(seur[['ATAC']]), pwm = pwm_set, genome = 'hg38',
    use.counts = FALSE
)
motif_obj = CreateMotifObject(data = motif_mat, pwm = pwm_set)
seur = SetAssayData(
    seur, assay = 'ATAC', layer = 'motifs', new.data = motif_obj
)
seur = RunChromVAR(object = seur, genome = BSgenome.Hsapiens.UCSC.hg38)

tripod_seur = getObjectsForModelFit(object = seur, chr = paste0("chr", 1:22))

seur = filterSeuratObject(object = seur, tripod.object = tripod_seur)
seur = processSeuratObject(
    object = seur, dim.rna = 1:50, dim.atac = 2:50, verbose = FALSE
)

################################################################################
#   Form metacells
################################################################################

cluster_df = optimizeResolution(
    object = seur,
    graph.name = "wsnn",
    assay.name = "WNN",
    resolutions = seq(10, 35, 5),
    min.num = 20
)

#   The vignette doesn't provide an algorithmic way to select the best
#   resolution, but provides this recommendation:
#       "Our empirical approach is to select a resolution that gives 80 or more
#        metacells, the majority of which contain 20 or more single cells"
#   I essentially implement that logic algorithmically here
best_res = cluster_df |>
    as_tibble() |>
    filter(num_clusters >= 80, num_below / num_clusters < 0.02) |>
    arrange(num_clusters) |>
    slice_head(n = 1) |>
    pull(resolution)

message(
    sprintf(
        "Selected resolution = %d. Here's the whole metric table:", best_res
    )
)
cluster_df |>
    as_tibble() |>
    print()

seur = getClusters(
    object = seur, graph.name = "wsnn", algorithm = 3, resolution = best_res,
    verbose = FALSE
)
metacell_seur = getMetacellMatrices(
    object = seur, cluster.name = "seurat_clusters"
)

################################################################################
#   Save in one list
################################################################################

pre_list = list(
    tripod_seur = tripod_seur,
    seur = seur,
    metacell metacell
)

qs_save(pre_list, out_path)

session_info()
