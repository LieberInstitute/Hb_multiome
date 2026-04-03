#   After computing linked peaks, there are a huge number of individual CSV
#   files for each combination of cell types, and there are duplicate results
#   (rows) across these files. Deduplicate, compute FDR from the deduplicated
#   p-values, and retain info about if the peak was called in the cell type
#   where the link was measured. Save a combined parquet file with this info,
#   as well as a version filtered by correlation and FDR thresholds

library(here)
library(tidyverse)
library(sessioninfo)
library(duckplyr)

cell_types = c(
    'Astrocyte', 'Endo', 'Excit.Thal', 'Inhib.Thal', 'LHb.1', 'LHb.1.3',
    'LHb.1.3.4', 'LHb.2.7', 'LHb.4', 'LHb.7', 'MHb.1', 'MHb.1.2', 'MHb.2',
    'MHb.3', 'Microglia', 'Oligo', 'OPC', 'Thal'
)
result_paths = here(
    'processed-data', '12_new_peaks', '01_link_peaks', '%s_%s.csv.gz'
)
full_out_path = here(
    'processed-data', '12_new_peaks', '01_link_peaks', 'all_data.parquet'
)
filtered_out_path = here(
    'processed-data', '12_new_peaks', '01_link_peaks', 'filtered_data.parquet'
)
cor_thres = 0.3
FDR_thres = 0.1

set.seed(0)
num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

#   First gather the key columns from every cell-type pair combination
result_df_list = list()
for (target_cell_type in cell_types) {
    for (other_cell_type in cell_types) {
        result_df_list[[length(result_df_list) + 1]] = read_csv_duckdb(
                sprintf(result_paths, target_cell_type, other_cell_type),
                prudence = "lavish"
            ) |>
            select(peak, gene, score, pvalue, target_cell_type, other_cell_type)
    }
}
result_df = bind_rows(result_df_list) |>
    #   Was the peak called in 'other_cell_type'?
    group_by(peak, gene) |>
    mutate(peak_called = other_cell_type %in% target_cell_type) |>
    ungroup() |>
    #   In many cases a link may have multiple rows indicating the peak was
    #   called in multiple cell types ('target_cell_type'), but we only care
    #   about each cell type where the link was measured ('other_cell_type') 
    distinct(peak, gene, other_cell_type, .keep_all = TRUE) |>
    dplyr::rename(cell_type = other_cell_type) |>
    group_by(cell_type) |>
    mutate(FDR = p.adjust(pvalue, method = "fdr")) |>
    ungroup() |>
    select(peak, gene, cell_type, score, FDR, peak_called) |>
    #   A combined version of the existing individual files, not filtered by
    #   correlation or FDR
    compute_parquet(full_out_path) |>
    #   Significance + correlation threshold
    filter(score > abs(cor_thres), FDR < FDR_thres) |>
    #   Is the link measured in multiple cell types?
    group_by(peak, gene) |>
    mutate(is_shared = length(unique(cell_type)) > 1) |>
    ungroup() |> 
    compute_parquet(filtered_out_path)

session_info()
