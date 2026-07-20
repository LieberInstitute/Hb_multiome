library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(qs2)
library(duckplyr)

seur_pb_path = here(
    'processed-data', '15_DARs', '01_pseudobulk_atac',
    'seur_pb_fine.qs2'
)
cell_map_path = here('raw-data', 'cell_type_map.csv')
dar_dir = here('processed-data', '15_DARs', '02_calculate_DARs')
out_dir = here("processed-data", "16_shiny_app", "04_prep_DARs")
log_fc_cutoff = log2(1.5)
p_adj_cutoff = 0.05

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(out_dir, showWarnings = FALSE)

cluster_map = read_csv(cell_map_path, show_col_types = FALSE)
rename_map <- stats::setNames(
    cluster_map$new_cell_type, cluster_map$old_cell_type
)

################################################################################
#   Prep DARs data frame
################################################################################

in_files = list.files(
    dar_dir, pattern = "^DARs_fine_.*\\.parquet$", full.names = TRUE
)
dar_df_list = list()
for (in_file in in_files) {
    dar_df_list[[in_file]] = read_parquet_duckdb(
            in_file, prudence = "lavish"
        ) |>
        select(peak, cell_type, avg_log2FC, p_val_adj) |>
        filter(p_val_adj < p_adj_cutoff, abs(avg_log2FC) > log_fc_cutoff) |>
        dplyr::rename(p_adj = p_val_adj)
}
dar_df = bind_rows(dar_df_list) |>
    collect() |>
    filter(str_detect(peak, '^chr')) |>
    mutate(
        DA_direction = factor(
            ifelse(avg_log2FC > 0, 'up', 'down'),
            levels = c('up', 'down')
        ),
        cell_type = dplyr::coalesce(
            unname(rename_map[cell_type]), cell_type
        )
    ) |>
    group_by(peak, DA_direction) |>
    arrange(p_adj) |>
    summarize(
        cell_type = paste(cell_type, collapse = ","),
        avg_log2FC = sign(avg_log2FC[1]) * max(abs(avg_log2FC)),
        p_adj = min(p_adj)
    ) |>
    ungroup()

dar_df = dar_df |>
    mutate(
        cell_type = factor(
            cell_type,
            levels = c(
                cluster_map$new_cell_type,
                setdiff(unique(dar_df$cell_type), cluster_map$new_cell_type)
            )
        )
    ) |>
    arrange(cell_type, DA_direction, p_adj, avg_log2FC) |>
    select(peak, cell_type, DA_direction, avg_log2FC, p_adj)

write_csv(dar_df, file.path(out_dir, "DARs.csv.gz"))

################################################################################
#   Prep psuedobulked ATAC data
################################################################################

seur_pb = qs_read(seur_pb_path)

#   Subset to 'data' layer and DAR peaks
seur_pb = DietSeurat(
    seur_pb, assays = "ATAC", layers = "data", misc = FALSE,
    features = unique(dar_df$peak)
)

seur_pb@meta.data$orig.ident = seur_pb[[]] |>
    as_tibble() |>
    mutate(
        orig.ident = str_replace_all(orig.ident, '-', '_'),
        orig.ident = dplyr::coalesce(
                unname(rename_map[orig.ident]), orig.ident
            ) |>
            factor(levels = cluster_map$new_cell_type)
    ) |>
    pull(orig.ident)
seur_pb@meta.data$mid_cluster = seur_pb@meta.data$orig.ident

qs_save(seur_pb, file = file.path(out_dir, "seur_pb_DARs.qs2"))

session_info()
