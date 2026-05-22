library(sessioninfo)
library(tidyverse)
library(here)
library(qs2)
library(Seurat)
library(Signac)
library(duckplyr)
library(cowplot)

prep_path = here(
    "processed-data", "13_tripod_trios", "02_tripod_preprocess",
    "preprocessed_objects_%s.qs2"
)
trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios.parquet"
)
num_examples = 5

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

trio_df = read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    filter(stringency_level == 1, is_unique, is_top_TF) |>
    select(-c(is_unique, stringency_level, is_top_TF)) |>
    collect()

#-------------------------------------------------------------------------------
#   Scatter plots: expression vs. accessibility, colored by TF expression
#-------------------------------------------------------------------------------

set.seed(8421)

sampled_trios = trio_df |>
    group_by(cell_type) |>
    slice_sample(n = num_examples) |>
    ungroup()

plot_list = list()
for (ct in unique(sampled_trios$cell_type)) {
    mats     = qs_read(sprintf(prep_path, ct))$metacell_seur
    ct_trios = sampled_trios |> filter(cell_type == ct)

    for (i in seq_len(nrow(ct_trios))) {
        row    = ct_trios[i, ]
        expr   = mats$rna[, row$gene, drop = TRUE]
        access = mats$peak[, row$peak, drop = TRUE]
        tf_expr = mats$rna[, row$TF, drop = TRUE]

        plot_list[[length(plot_list) + 1]] = ggplot(
                data.frame(access = access, expr = expr, tf_expr = tf_expr),
                aes(x = access, y = expr, color = tf_expr)
            ) +
            geom_point(size = 0.15) +
            scale_color_viridis_c(name = "TF expr.") +
            labs(
                title = sprintf("%s\n%s | %s", ct, row$gene, row$TF),
                x = "Access.", y = "Expr."
            ) +
            theme_bw(base_size = 3.5) +
            theme(
                plot.title  = element_text(size = 2.5, lineheight = 1.1),
                legend.position = "none"
            )
    }

    rm(mats); gc()
}

plot_grid(plotlist = plot_list, ncol = num_examples)

session_info()
