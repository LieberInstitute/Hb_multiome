#   Generate one CSV of file paths and MD5 sums for processed data

library(tidyverse)
library(here)
library(sessioninfo)
library(spatialLIBD)

visium_repo_dir = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium'
multiome_map_path = here('raw-data', 'sample_id_map.csv')
out_path = here(
    'processed-data', '19_data_uploads', '02_processed_data', 'map.csv'
)
out_dir = here('processed-data', '19_data_uploads', '01_file_map', 'fastq_flat')
fetch_object_names = c(
    "habenula_atlas_HD_spe_cell", "habenula_atlas_HD_spe_cell_pseudobulk",
    "habenula_atlas_visium_spe", "habenula_atlas_visium_spe_pseudobulk",
    "habenula_atlas_snMultiome_seurat_cell",
    "habenula_atlas_snMultiome_seurat_metacell"
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)

################################################################################
#   RDS files
################################################################################

rds_paths = c()
for (object_name in fetch_object_names) {
    x = fetch_data(object_name)
    file_path = here(out_dir, paste0(object_name, '.rds'))
    saveRDS(x, file = file_path)
    rds_paths = c(rds_paths, file_path)
}

################################################################################
#   Spaceranger processed data
################################################################################

#-------------------------------------------------------------------------------
#   Visium
#-------------------------------------------------------------------------------

session_info()
  