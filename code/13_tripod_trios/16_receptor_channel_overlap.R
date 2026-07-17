library(here)
library(tidyverse)
library(sessioninfo)
library(duckplyr)

anno_path = here(
    'processed-data', '05_03_annotation_adjustments', '16_ion_and_receptor_exp',
    'genes_of_interest.csv'
)
trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

overlap_df = read_parquet_duckdb(trio_path, prudence = "stingy") |>
    dplyr::filter(!is_intersect | (is_intersect & stringency_level == 1)) |>
    collect() |>
    mutate(
        stringency_level = ifelse(is_intersect, 'Intersect', stringency_level)
    ) |>
    dplyr::select(peak, gene, TF, cell_type, stringency_level) |>
    inner_join(read_csv(anno_path, show_col_types = FALSE), by = 'gene')

overlap_df |>
    filter(stringency_level == 'Intersect', grepl('[ML]Hb', cell_type))
