#   After computing linked peaks, there are a large number of individual
#   parquet files with both extra and missing information. Merge all files,
#   drop unused columns, and add info about if the peak was called in the
#   given cell type. Save an unfiltered version and version filtered by
#   correlation and FDR

library(here)
library(tidyverse)
library(sessioninfo)
library(duckplyr)

dataset = c("pb", "metacell")[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

if (dataset == "pb") {
    result_paths = here(
        'processed-data', '12_new_peaks', '01_link_peaks', '%s.parquet'
    )
} else {
    result_paths = here(
        "processed-data", "12_new_peaks", "14_metacell_link_peaks", "%s.parquet"
    )
}

cell_types = c(
    "Astrocyte", "Endo", "Ependymal", "Excit.Thal", "Inhib_LHb_4.1",
    "Inhib_LHb_4.2", "Inhib.Thal", "LHb.1.3.4", "LHb.2.7", "LHb.4", "MHb.1",
    "MHb.1.2", "MHb.2", "MHb.3", "Microglia", "Oligo", "OPC"
)
peak_path = here(
    'processed-data', '11_link_prep', '01_call_peaks',
    'macs3_peaks.csv.gz'
)
full_out_path = here(
    'processed-data', '12_new_peaks', '17_filter_links',
    sprintf('%s_all_data.parquet', dataset)
)
filtered_out_path = here(
    'processed-data', '12_new_peaks', '17_filter_links',
    sprintf('%s_filtered_data.parquet', dataset)
)
cor_thres = 0.3

set.seed(0)
num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(full_out_path), showWarnings = FALSE)

peak_df = read_csv_duckdb(peak_path, prudence = "lavish") |>
    #   Half to be careful about formatting here because I encountered a rare
    #   case where one peak end value, 200000, was converted to 2e+05 in the
    #   default paste call (before using formatC later), leading to failed joins
    mutate(
        peak = paste(
            seqnames,
            formatC(start, format = "d"),
            formatC(end, format = "d"),
            sep = "-"
        )
    ) |>
    select(peak, peak_called_in) |>
    collect()

#   First gather the key columns from every cell type
result_df_list = list()
for (this_cell_type in cell_types) {
    result_df_list[[length(result_df_list) + 1]] = read_parquet_duckdb(
            sprintf(result_paths, this_cell_type), prudence = "lavish"
        ) |>
        select(peak, gene, cell_type, score, FDR)
}
result_df = bind_rows(result_df_list) |>
    left_join(peak_df, by = "peak") |>
    mutate(
        peak_called = str_detect(
            peak_called_in, sprintf('(^|,)%s(,|$)', cell_type)
        )
    ) |>
    select(peak, gene, cell_type, score, FDR, peak_called) |>
    #   A combined version of the existing individual files, not filtered by
    #   correlation or FDR
    compute_parquet(full_out_path) |>
    #   Correlation threshold only
    filter(abs(score) > cor_thres) |>
    #   Is the link measured in multiple cell types?
    group_by(peak, gene) |>
    mutate(is_shared = length(unique(cell_type)) > 1) |>
    ungroup() |> 
    compute_parquet(filtered_out_path)

session_info()
