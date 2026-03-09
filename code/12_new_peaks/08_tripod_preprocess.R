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

seur_path = here(
    "processed-data", "12_new_peaks", "06_non_pb_seur",
    "non_pb_seur.qs2"
)

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
tx_gr = tripod.pbmc$transcripts.gr
peaks_gr = tripod.pbmc$peaks.gr
motif_x_tf = tripod.pbmc$motifxTF
peak_x_motif = tripod.pbmc$peakxmotif

seur = filterSeuratObject(object = seur, tripod.object = tripod_seur)