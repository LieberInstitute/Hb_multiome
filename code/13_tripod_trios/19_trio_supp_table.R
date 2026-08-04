library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
out_path = here(
    "processed-data", "13_tripod_trios", "19_trio_supp_table",
    "trios.csv"
)
cell_map_path = here('raw-data', 'cell_type_map.csv')

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(out_path), showWarnings = FALSE)

cluster_map = read_csv(cell_map_path, show_col_types = FALSE)
rename_map <- stats::setNames(
    cluster_map$new_cell_type, cluster_map$old_cell_type
)

read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    dplyr::rename(p_adj = adj) |>
    select(peak, gene, TF, coef, p_adj, stringency_level, cell_type) |>
    collect() |>
    pivot_wider(
        names_from = stringency_level,
        values_from = c(coef, p_adj),
        names_glue = "{.value}_level{stringency_level}"
    ) |>
    mutate(
        cell_type = unname(rename_map[cell_type]),
        test_level = dplyr::case_when(
            !is.na(coef_level1) & !is.na(coef_level2) ~ 'Intersect',
            !is.na(coef_level1) ~ '1',
            !is.na(coef_level2) ~ '2',
            TRUE ~ NA
        )
    ) |>
    select(
        cell_type, test_level, peak, gene, TF, coef_level1, p_adj_level1,
        coef_level2, p_adj_level2
    ) |>
    arrange(test_level, cell_type, desc(coalesce(coef_level1, coef_level2))) |>
    write_csv(out_path)

session_info()
