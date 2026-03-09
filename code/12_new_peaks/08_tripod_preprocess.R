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

seur = qs_read(seur_path)

DefaultAssay(seur) = "ATAC"
pwm_set = getMatrixSet(
    x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE)
)
motif_mat = CreateMotifMatrix(
    features = granges(seur), pwm = pwm_set, genome = 'hg38', use.counts = FALSE
)
motif_obj = CreateMotifObject(data = motif_mat, pwm = pwm_set)
seur = SetAssayData(seur, assay = 'ATAC', layer = 'motifs', new.data = motif_obj)
seur = RunChromVAR(object = seur, genome = BSgenome.Hsapiens.UCSC.hg38)

tripod.pbmc <- getObjectsForModelFit(object = pbmc, chr = paste0("chr", 1:22))
transcripts.gr <- tripod.pbmc$transcripts.gr
peaks.gr <- tripod.pbmc$peaks.gr
motifxTF <- tripod.pbmc$motifxTF
peakxmotif <- tripod.pbmc$peakxmotif
