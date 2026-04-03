#   Call peaks, providing cell types as the grouping variable
#
#   It appears that MACS3 is usable with CallPeaks from Signac 1.17.0
#   (see https://github.com/stuart-lab/signac/issues/1085#issuecomment-2402690609).
#   Note that I didn't do the officially recommended solution
#   (https://github.com/stuart-lab/signac/issues/1085#issuecomment-4115198769),
#   which uses the v2 (beta) branch of the github version of Signac, as this
#   version introduces breaking changes that prevent interacting with the Seurat
#   object imported in this script. I figured it would be better to use the stable
#   version, with the hack of providing the MACS3 path to the macs2.path parameter
#   of CallPeaks. To install MACS3, I ran:
#       uv pip install --upgrade macs3
#   The resulting binary, despite being installed in a python environment,
#   appears to work as a standalone tool (maybe it could be made into a module?)

library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(qs2)

macs3_path = '/users/neagles/.conda/envs/macs3_env/bin/macs3'
seur_path = here(
    'processed-data', '05_03_annotation_adjustments', '06_refined_annotations',
    'refined_annotation_multiomeHab_Seurat.qs2'
)
out_path = here(
    'processed-data', '11_link_prep', '01_call_peaks',
    'macs3_peaks.csv.gz'
)

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

seur = qs_read(seur_path)
DefaultAssay(seur) = "ATAC"

peaks = CallPeaks(
    object = seur, group.by = 'refined_mid_cluster', macs2.path = macs3_path,
    verbose = TRUE
)

write_csv(peaks, out_path)

session_info()
