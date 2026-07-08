library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
dar_dir = here("processed-data", "15_DARs", "02_calculate_DARs")
out_path = here(
    "processed-data", "15_DARs", "04_trio_overlap",
    "trio_dar_overlap.csv"
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

trio_df = read_parquet_duckdb(trio_path, prudence = "stingy") |>
    filter(is_intersect, stringency_level == 1) |>
    select(peak, gene, TF, cell_type, coef, adj) |>
    collect()

dar_df = lapply(
        list.files(
            dar_dir, pattern = '^DARs_(fine|mid)_[ML]Hb.*\\.parquet$',
            full.names = TRUE
        ),
        function(x) read_parquet_duckdb(x, prudence = 'lavish')
    ) |>
    bind_rows() |>
    filter(p_val_adj < 0.1) |>
    select(peak, p_val_adj, avg_log2FC, cell_type, resolution) |>
    dplyr::rename(
        dar_p_adj = p_val_adj, dar_log_fc = avg_log2FC,
        dar_cell_type = cell_type
    ) |>
    collect() |>
    mutate(broad_cell_type = str_extract(dar_cell_type, "^[ML]Hb"))

trio_df |>
    filter(str_detect(cell_type, "^[ML]Hb")) |>
    mutate(broad_cell_type = str_extract(cell_type, "^[ML]Hb")) |>
    dplyr::rename(
        trio_cell_type = cell_type, trio_p_adj = adj, trio_coef = coef
    ) |>
    inner_join(dar_df, by = c("peak", "broad_cell_type")) |>
    select(
        peak, gene, TF, trio_cell_type, dar_cell_type, trio_coef, trio_p_adj,
        dar_p_adj, dar_log_fc
    ) |>
    write_csv(out_path)
  
session_info()
