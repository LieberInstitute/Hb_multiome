#   Generate one CSV of file paths and MD5 sums for processed data

library(tidyverse)
library(here)
library(sessioninfo)
library(spatialLIBD)

visium_repo_dir = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium'
multiome_map_path = here('raw-data', 'sample_id_map.csv')
cellranger_dir = here('processed-data', 'cellrangerARC')
cellranger_reanalyze_dir = here('processed-data', 'cellrangerARC_reanalyze')
hd_spaceranger_dir = here(
    visium_repo_dir, 'processed-data', '01_spaceranger', 'five_samples_10_2025'
)
he_spaceranger_dir = here(visium_repo_dir, 'processed-data', '01_spaceranger')
out_path = here(
    'processed-data', '19_data_uploads', '02_processed_data', 'map.csv'
)
out_dir = here('processed-data', '19_data_uploads', '01_file_map', 'fastq_flat')
rds_dir = here('processed-data', '19_data_uploads', '02_processed_data', 'rds')
fetch_object_names = c(
    "habenula_atlas_HD_spe_cell", "habenula_atlas_HD_spe_cell_pseudobulk",
    "habenula_atlas_visium_spe", "habenula_atlas_visium_spe_pseudobulk",
    "habenula_atlas_snMultiome_seurat_cell",
    "habenula_atlas_snMultiome_seurat_metacell"
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)
dir.create(rds_dir, showWarnings = FALSE)

################################################################################
#   RDS files
################################################################################

rds_processed_df = fetch_object_names |>
    map(
        function(object_name) {
            x = fetch_data(object_name)
            file_path = here(rds_dir, paste0(object_name, '.rds'))
            saveRDS(x, file = file_path)
            tibble(object_name = object_name, file_path = file_path)
        }
    ) |>
    list_rbind() |>
    mutate(file_name = basename(file_path))

################################################################################
#   Spaceranger processed data
################################################################################

#-------------------------------------------------------------------------------
#   Visium HD
#-------------------------------------------------------------------------------

hd_sample_dirs = list.dirs(hd_spaceranger_dir, full.names = FALSE, recursive = FALSE)

hd_subdirs = c('raw_feature_bc_matrix', 'spatial')

hd_processed_paths = hd_sample_dirs |>
    map(
        function(sample_dir) {
            donor = sprintf(
                'Br%s', str_extract(sample_dir, '[AD]1_([0-9]{4})', group = 1)
            )

            hd_subdirs |>
                map(
                    function(subdir) {
                        dir_path = here(
                            hd_spaceranger_dir, sample_dir, 'outs',
                            'binned_outputs', 'square_002um', subdir
                        )
                        file_path = list.files(
                            dir_path, full.names = TRUE, recursive = TRUE
                        )
                        tibble(donor = donor, file_path = file_path)
                    }
                ) |>
                list_rbind()
        }
    ) |>
    list_rbind()

hd_processed_df = hd_processed_paths |>
    #   This file is redundant with 'tissue_positions.parquet' and isn't
    #   needed for upload
    filter(basename(file_path) != 'tissue_positions.csv') |>
    mutate(
        file_path = normalizePath(file_path),
        file_name = paste(donor, basename(file_path), sep = '_')
    )

#-------------------------------------------------------------------------------
#   Visium (standard)
#-------------------------------------------------------------------------------

#   Match the exact set of capture areas whose raw data is being uploaded
spe = fetch_data('habenula_atlas_visium_spe')
he_capture_areas = unique(spe$sample_id)

he_subdirs = c('raw_feature_bc_matrix', 'spatial')

he_processed_paths = he_capture_areas |>
    map(
        function(capture_area) {
            he_subdirs |>
                map(
                    function(subdir) {
                        dir_path = here(
                            he_spaceranger_dir, capture_area, 'outs', subdir
                        )
                        file_path = list.files(
                            dir_path, full.names = TRUE, recursive = TRUE
                        )
                        tibble(capture_area = capture_area, file_path = file_path)
                    }
                ) |>
                list_rbind()
        }
    ) |>
    list_rbind()

he_processed_df = he_processed_paths |>
    #   Unlike Visium HD, there's no 'tissue_positions.parquet' alternative
    #   here, so 'tissue_positions.csv' is kept
    mutate(
        file_path = normalizePath(file_path),
        #   Capture areas (not donors) are the unit that needs
        #   deduplication: each donor has 4 capture areas, so donor alone
        #   isn't a unique prefix
        file_name = paste(capture_area, basename(file_path), sep = '_')
    )

################################################################################
#   Cellranger-ARC processed data (multiome)
################################################################################

#   'sample_id_2' is the Cell Ranger ARC sample-directory naming convention
#   (distinct from 'sample_id_1', used for raw FASTQs) and is the exact set
#   of samples whose raw data is being uploaded
multiome_map_df = read_csv(multiome_map_path, show_col_types = FALSE)

cellranger_processed_paths = multiome_map_df$sample_id_2 |>
    map(
        function(cr_sample) {
            #   Most QC-relevant outputs come from the original (non-reanalyzed)
            #   Cell Ranger ARC run
            cr_dir = here(cellranger_dir, cr_sample, 'outs')
            cr_files = file.path(
                cr_dir, c('atac_fragments.tsv.gz', 'per_barcode_metrics.csv')
            )

            #   The filtered feature-barcode matrix comes from the reanalyzed
            #   run instead, which applies an updated cell-calling algorithm
            reanalyze_dir = here(
                cellranger_reanalyze_dir, paste0(cr_sample, '_reanalysis'), 'outs'
            )
            reanalyze_files = file.path(
                reanalyze_dir, 'filtered_feature_bc_matrix.h5'
            )

            tibble(sample_id_2 = cr_sample, file_path = c(cr_files, reanalyze_files))
        }
    ) |>
    list_rbind()

stopifnot(
    'some expected Cell Ranger ARC processed-data files are missing' =
        all(file.exists(cellranger_processed_paths$file_path))
)

cellranger_processed_df = cellranger_processed_paths |>
    #   Match the canonical 'sample_id' naming convention established in
    #   01_file_map.R (e.g. 'S03_Hb', not the raw directory name 'S3_Hb_KDM'
    #   or '4S_Hb_KDM') so symlink names are consistent across the upload
    left_join(multiome_map_df, by = 'sample_id_2') |>
    mutate(
        sample_id = sub('_r$', '', sample_id_1),
        file_path = normalizePath(file_path),
        #   Deduplicate by sample, since e.g. 'per_barcode_metrics.csv' is
        #   shared across every sample's directory structure
        file_name = paste(sample_id, basename(file_path), sep = '_')
    ) |>
    select(-sample_id_1, -sample_id_2)

################################################################################
#   Combine, symlink, compute checksums, and write out
################################################################################

processed_df = bind_rows(
    rds_processed_df, hd_processed_df, he_processed_df, cellranger_processed_df
)

stopifnot(
    'file_name is not unique across the combined processed data' =
        !any(duplicated(processed_df$file_name))
)

flat_path = file.path(out_dir, processed_df$file_name)
unlink(flat_path)
file.symlink(processed_df$file_path, flat_path)

processed_df |>
    mutate(file_path = flat_path) |>
    mutate(md5_checksum = tools::md5sum(file_path)) |>
    select(file_path, md5_checksum) |>
    write_csv(out_path)

session_info()
  