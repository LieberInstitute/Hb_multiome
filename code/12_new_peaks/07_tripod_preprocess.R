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

seur_path = here(
    "processed-data", "12_new_peaks", "06_non_pb_seur",
    "non_pb_seur.qs2"
)

seur = qs_read(seur_path)
