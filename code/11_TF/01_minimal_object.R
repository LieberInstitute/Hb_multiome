#   It's taking too long to interactively test TF code because the full Seurat
#   object is massive. Drop unnecessary assays and reduced dimensions and save
#   a minimal object

library(Seurat)
library(Signac)
library(here)
library(lobstr)
library(sessioninfo)

in_path = here(
    'processed-data', '06_peak_calling', '11_peaks_merge_MACS2',
    'Seurat_peaks_merged_cell_level_Mid_resolution.rds'
)
out_path = here('processed-data', '11_TF', 'minimal_seur.rds')

seur = readRDS(in_path)
message("Original object size: ", format(obj_size(seur), units = "GB"))

#   Drop assays
seur[['RNA']] = NULL
seur[['ATAC']] = NULL
DefaultAssay(seur) = "ATAC_macs2_merged"

#   Drop reduced dims
seur@reductions = list()

message("Reduced object size: ", format(obj_size(seur), units = "GB"))

saveRDS(seur, out_path)

session_info()
