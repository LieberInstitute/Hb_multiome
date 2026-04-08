library(sessioninfo)
library(tidyverse)
library(here)
library(ComplexHeatmap)
library(circlize)
library(duckplyr)
library(Seurat)
library(Signac)
library(qs2)

trio_paths = here(
    "processed-data", "13_tripod_trios", "03_tripod_trios", "trios_%s.parquet"
)
link_path = here(
    "processed-data", "12_new_peaks", "01_link_peaks", "filtered_data.parquet"
)
seur_path = here(
    'processed-data', '11_link_prep', '02_rebuild_atac_assay',
    'cell_level_seur.qs2'
)
plot_path = here(
    "plots", "13_tripod_trios", "04_trio_heatmap", "heatmap.pdf"
)
cell_type1 = "MHb.2"
cell_type2 = "LHb.2.7"
cor_thres = 0.1
FDR_thres = 0.1

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(plot_path), showWarnings = FALSE)

link_df = read_parquet_duckdb(link_path, prudence = 'stingy') |>
    collect()

trio_df_list = list()
for (cell_type in c(cell_type1, cell_type2)) {
    trio_df_list[[cell_type]] = read_parquet_duckdb(
            sprintf(trio_paths, cell_type), prudence = 'lavish'
        ) |>
        select(peak, gene, TF, adj, condition_on, stringency_level) |>
        mutate(cell_type = cell_type)
}
trio_df = bind_rows(trio_df_list) |>
    collect()

final_df = trio_df |>
    filter(condition_on == 'Yj', stringency_level == 2) |>
    select(peak, gene, TF, adj, cell_type) |>
    inner_join(
        link_df |>
            dplyr::rename(cell_type = target_cell_type) |>
            select(peak, gene, cell_type),
        by = c("peak", "gene", "cell_type")
    )

seur = qs_read(seur_path)

CoveragePlot(seur, region = final_df$peak[1], features = final_df$gene[1])
