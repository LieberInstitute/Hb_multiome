library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(qs2)
library(GenomicRanges)
library(duckplyr)

#   Grab this specific cell type and resolution from the array task
task_id = as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

task_map_path = here(
    'processed-data', '15_DARs', '01_pseudobulk_atac',
    'task_map.csv'
)
task_map = read_csv(task_map_path, show_col_types = FALSE)
resolution = task_map$resolution[task_id]
cell_type = task_map$cell_type[task_id]

seur_path = here(
    'processed-data', '15_DARs', '01_pseudobulk_atac',
    sprintf('seur_pb_%s.qs2', resolution)
)
out_path = here(
    'processed-data', '15_DARs', '02_calculate_DARs',
    sprintf('DARs_%s_%s.parquet', resolution, cell_type)
)
p_adj_cutoff = 0.2

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(out_path), showWarnings = FALSE)

seur_pb = qs_read(seur_path)

#   Fix cell types
seur_pb@meta.data$orig.ident = str_replace_all(
    seur_pb@meta.data$orig.ident, "-", "_"
)
Idents(seur_pb) = seur_pb$orig.ident

#   1-vs-all approach with minimal filtering. Roughly based off of:
#   https://stuartlab.org/signac/articles/pbmc_vignette.html#find-differentially-accessible-peaks-between-cell-types
temp = FindMarkers(
        object = seur_pb, ident.1 = cell_type, ident.2 = NULL, only.pos = FALSE,
        logfc.threshold = 0, min.pct = 0.1
    ) |>
    rownames_to_column("peak") |>
    as_tibble() |>
    filter(p_val_adj < p_adj_cutoff) |>
    mutate(cell_type = cell_type, resolution = resolution) |>
    compute_parquet(out_path)

session_info()
