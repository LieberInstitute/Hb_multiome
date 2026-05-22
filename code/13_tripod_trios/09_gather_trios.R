library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)

cell_types = c(
    "Astrocyte", "Endo", "Ependymal", "Excit.Thal", "Inhib_LHb_4.1",
    "Inhib_LHb_4.2", "Inhib.Thal", "LHb.1.3.4", "LHb.2.7", "LHb.4", "MHb.1",
    "MHb.1.2", "MHb.2", "MHb.3", "Microglia", "Oligo", "OPC",
    "MHb", "LHb", "Inhib_LHb"
)
trio_paths = here(
    "processed-data", "13_tripod_trios", "03_tripod_trios", "trios_%s.parquet"
)
out_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios.parquet"
)
fdr_thres = 0.05
background_fdr_thres = 0.2

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

# Read trios at FDR < 0.05 (main set) and FDR < 0.2 (background for uniqueness)
trio_df_list = list()
background_df_list = list()
for (cell_type in cell_types) {
    this_path = sprintf(trio_paths, cell_type)
    if (file.exists(this_path)) {
        raw = read_parquet_duckdb(this_path, prudence = 'lavish') |>
            filter(condition_on == 'Yj', adj < background_fdr_thres) |>
            select(peak, gene, TF, coef, adj, stringency_level) |>
            mutate(cell_type = cell_type) |>
            collect()
        background_df_list[[cell_type]] = raw
        trio_df_list[[cell_type]] = raw |> filter(adj < fdr_thres)
    } else {
        message(
            sprintf(
                "Trios were not found for cell type %s; skipping", cell_type
            )
        )
    }
}

# Build set of peak-gene pairs per cell type from the background
background_pg_list = lapply(background_df_list, \(df) {
    unique(paste(df$peak, df$gene, sep = "|"))
})

# Export final tibble at FDR < 0.05, annotated with is_unique and is_intersect
bind_rows(trio_df_list) |>
    mutate(pg = paste(peak, gene, sep = "|")) |>
    group_by(cell_type) |>
    mutate(
        # is_unique: peak-gene not seen in any other cell type at FDR < 0.2
        is_unique = !(
            pg %in% unlist(
                background_pg_list[names(background_pg_list) != cell_type[1]]
            )
        ),
        # is_intersect: within cell type, trio (by peak, gene, TF) appears at
        # both stringency levels
        pg_tf = paste(peak, gene, TF, sep = "|"),
        is_intersect = pg_tf %in% pg_tf[stringency_level == 1] &
            pg_tf %in% pg_tf[stringency_level == 2]
    ) |>
    #   Finally annotate the most significant TF for peak-gene pairs with
    #   multiple TFs
    group_by(peak, gene, cell_type, stringency_level) |>
    arrange(adj) |>
    mutate(is_top_TF = row_number() == 1) |>
    ungroup() |>
    #   Clean up and export
    arrange(cell_type, stringency_level, adj) |>
    select(-pg, -pg_tf) |>
    compute_parquet(out_path)

session_info()
