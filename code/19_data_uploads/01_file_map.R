#   Generate one CSV mapping between donor, sample ID, library ID, and file path
#   for all files to upload

library(tidyverse)
library(here)
library(sessioninfo)

multiome_map_path = here('raw-data', 'sample_id_map.csv')
multiome_fastq_dir1 = here('raw-data', 'FASTQ_2024')
multiome_fastq_dir2 = here('raw-data', 'FASTQ_2024_data_package2')

multiome_map_df = read_csv(multiome_map_path, show_col_types = FALSE)

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
        file_path = normalizePath(file_path)
    ) |>
    select(donor, sample_id, library_id, file_path)


