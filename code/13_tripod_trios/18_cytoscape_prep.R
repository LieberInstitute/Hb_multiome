library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
out_path = here(
    "processed-data", "13_tripod_trios", "18_cytoscape_prep",
    "trio_cytoscape.csv"
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(out_path), showWarnings = FALSE)

read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    filter(is_intersect, stringency_level == 1) |>
    select(gene, TF, cell_type, coef, adj) |>
    collect() |>
    write_csv(out_path)

session_info()
