#   Generate one CSV mapping between donor, sample ID, library ID, and file path
#   for all files to upload

library(tidyverse)
library(here)
library(sessioninfo)

visium_repo_dir = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium'
multiome_map_path = here('raw-data', 'sample_id_map.csv')
multiome_fastq_dir1 = here('raw-data', 'FASTQ_2024')
multiome_fastq_dir2 = here('raw-data', 'FASTQ_2024_data_package2')
hd_fastq_dir = here(visium_repo_dir, 'raw-data', 'fastqs')
hd_image_dir = here(visium_repo_dir, 'raw-data', 'images', 'vis-hd')
hd_image_info_path = here(
    visium_repo_dir, 'code', '01_spaceranger', 'all_hd_samples_10_2025.txt'
)

demo_path = here(
    visium_repo_dir, 'processed-data', '14_supp_tables', 'donor_demographics.csv'
)

################################################################################
#   Multiome data
################################################################################

multiome_map_df = read_csv(multiome_map_path, show_col_types = FALSE)

#-------------------------------------------------------------------------------
#   FASTQs
#-------------------------------------------------------------------------------

multiome_fastq = c(
    list.files(
        multiome_fastq_dir1, pattern = 'fastq.gz$', full.names = TRUE,
        recursive = TRUE
    ),
    list.files(
        multiome_fastq_dir2, pattern = 'fastq.gz$', full.names = TRUE,
        recursive = TRUE
    )
)

multiome_fastq_df = tibble(file_path = multiome_fastq) |>
    mutate(
        sample_id = sprintf(
            'S%02d_Hb_r',
            as.integer(str_extract(basename(file_path), '^[0-9]+'))
        ),
        library_id = sprintf(
            'lib_%s_%s',
            str_extract(sample_id, 'S[0-9]{2}_Hb'),
            str_extract(file_path, 'GEX|ATAC')
        )
    ) |>
    left_join(multiome_map_df, by = c('sample_id' = 'sample_id_1')) |>
    mutate(
        sample_id = sub('_r$', '', sample_id),
        file_path = normalizePath(file_path),
        open_access = FALSE
    ) |>
    select(donor, sample_id, library_id, file_path, open_access)

################################################################################
#   Visium HD data
################################################################################

#-------------------------------------------------------------------------------
#   FASTQs
#-------------------------------------------------------------------------------

hd_fastq = list.files(
    hd_fastq_dir, pattern = 'fastq.gz$', full.names = TRUE, recursive = TRUE
)

hd_fastq_df = tibble(file_path = hd_fastq) |>
    filter(
        grepl('[AD]1_[0-9]{4}', file_path),
        !grepl('[AD]1_(9037|8518)', file_path)
    ) |>
    mutate(
        donor = sprintf(
            'Br%s', str_extract(file_path, '[AD]1_([0-9]{4})', group = 1)
        ),
        sample_id = donor,
        library_id = paste('lib', donor, sep = '_'),
        file_path = normalizePath(file_path),
        open_access = FALSE
    ) |>
    select(donor, sample_id, library_id, file_path, open_access)

#-------------------------------------------------------------------------------
#   Images
#-------------------------------------------------------------------------------

hd_image_df = read_table(
        hd_image_info_path, show_col_types = FALSE,
        col_names = c('sample_id', 'image_id')
    ) |>
    mutate(
        donor = sprintf('Br%s', str_extract(sample_id, '[0-9]{4}$')),
        sample_id = donor,
        library_id = paste('lib', donor, sep = '_'),
        file_path = file.path(hd_image_dir, paste0(image_id, '.tif')) |>
            normalizePath(),
        open_access = TRUE
    ) |>
    select(donor, sample_id, library_id, file_path, open_access)
