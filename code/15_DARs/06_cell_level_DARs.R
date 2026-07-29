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
cell_type = task_map$cell_type[task_id]

seur_path = here(
    'processed-data', '11_link_prep', '02_rebuild_atac_assay',
    'cell_level_seur.qs2'
)
out_path = here(
    'processed-data', '15_DARs', '06_cell_level_DARs',
    sprintf('DARs_%s.parquet', cell_type)
)
p_adj_cutoff = 0.2

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(out_path), showWarnings = FALSE)

seur = qs_read(seur_path)

Idents(seur) = seur$refined_mid_cluster

#   1-vs-all approach with minimal filtering. Roughly based off of:
#   https://stuartlab.org/signac/articles/pbmc_vignette.html#find-differentially-accessible-peaks-between-cell-types
temp = FindMarkers(
        object = seur, slot = "data", ident.1 = cell_type, ident.2 = NULL,
        logfc.threshold = 0
    ) |>
    rownames_to_column("peak") |>
    as_tibble() |>
    filter(p_val_adj < p_adj_cutoff) |>
    mutate(cell_type = cell_type) |>
    compute_parquet(out_path)

session_info()
