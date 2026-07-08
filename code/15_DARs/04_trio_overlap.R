library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
dar_path = here(
    "processed-data", "15_DARs", "03_gather", "DARs_fine.csv.gz"
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

trio_df = read_parquet_duckdb(trio_path, prudence = "stingy") |>
    filter(is_intersect, stringency_level == 1) |>
    select(peak, gene, TF, cell_type, coef, adj) |>
    collect()

dar_df = read_csv(dar_path, show_col_types = FALSE) |>
    dplyr::rename(dar_p_adj = p_val_adj, dar_log_fc = avg_log2FC) |>
    select(peak, cell_type, dar_p_adj, dar_log_fc)

message("Intersection trios with DA peak:")
trio_df |>
    inner_join(dar_df, by = c("peak", "cell_type")) |>
    print()

session_info()
