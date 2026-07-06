library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)

in_dir = here('processed-data', '15_DARs', '02_calculate_DARs')
out_dir = here('processed-data', '15_DARs', '03_gather')
p_adj_cutoff = 0.05
log_fc_cutoff = log2(3)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(out_dir, showWarnings = FALSE)

in_files = list.files(
    in_dir, pattern = "\\.parquet$", full.names = TRUE
)
dar_df_list = list()
for (in_file in in_files) {
    dar_df_list[[in_file]] = read_parquet_duckdb(
            in_file, prudence = "lavish"
        ) |>
        select(peak, cell_type, resolution, avg_log2FC, p_val_adj) |>
        filter(p_val_adj < p_adj_cutoff, avg_log2FC > log_fc_cutoff)
}
dar_df = bind_rows(dar_df_list) |>
    collect()

dar_df |>
    filter(!str_detect(resolution, '^[ML]Hb|Thal')) |>
    write_csv(file.path(out_dir, 'DARs_broad.csv.gz'))

dar_df |>
    filter(
        (resolution == 'mid') |
        ((resolution == 'fine') & !str_detect(cell_type, '^[ML]Hb|Thal'))
    ) |>
    write_csv(file.path(out_dir, 'DARs_mid.csv.gz'))

dar_df |>
    filter(resolution == 'fine') |>
    write_csv(file.path(out_dir, 'DARs_fine.csv.gz'))

session_info()
