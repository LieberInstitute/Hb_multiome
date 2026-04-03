#   At the time of creating this script, in order to use MACS3 rather than
#   MACS2, I had to:
#       in R: remotes::install_github('stuart-lab/signac@v2')
#       in terminal: uv pip install --upgrade macs3
#   The MACS3 binary seems to work without needing the specific Python
#   environment loaded, so maybe this could be made into a module to more
#   easily reproduce this code as a different user?
#
#   Call peaks, providing cell types as the grouping variable

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
    object = seur, group.by = 'refined_mid_cluster', macs3.path = macs3_path,
    verbose = TRUE
)

write_csv(peaks, out_path)

session_info()
