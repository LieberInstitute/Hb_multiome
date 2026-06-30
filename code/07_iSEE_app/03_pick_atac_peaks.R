library(here)
library(tidyverse)
library(qs2)
library(sessioninfo)
library(SingleCellExperiment)
library(duckplyr)

sce_path = here(
    'processed-data', '07_iSEE_app', '01_prep_sce', 'sce_ATAC_iSEE.qs2'
)
trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
peak_path = here(
    'processed-data', '11_link_prep', '01_call_peaks',
    'macs3_peaks.csv.gz'
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

peak_df = read_csv_duckdb(peak_path, prudence = "stingy") |>
    select(seqnames, start, end, peak_called_in) |>
    collect() |>
    mutate(peak = sprintf("%s-%d-%d", seqnames, start, end)) |>
    select(peak, peak_called_in)

trio_df = read_parquet_duckdb(trio_path, prudence = "stingy") |>
    filter(is_intersect, stringency_level == 1) |>
    collect() |>
    left_join(peak_df, by = 'peak')

stopifnot(all(trio_df$peak %in% rownames(sce)))

marker_df = findMarkers_1vAll(
        sce[unique(trio_df$peak),], assay_name = 'logcounts',
        cellType_col = 'mid_cluster'
    ) |>
    dplyr::rename(
        peak = gene, marker_fdr = log.FDR, cell_type = cellType.target
    ) |>
    select(peak, marker_fdr, cell_type)

trio_df |>
    left_join(marker_df, by = c('peak', 'cell_type')) |>
    group_by(cell_type) |>
    arrange(marker_fdr) |>
    slice_head(n = 1) |>
    select(peak, cell_type)
