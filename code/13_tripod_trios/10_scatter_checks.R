library(sessioninfo)
library(tidyverse)
library(here)
library(qs2)
library(Seurat)
library(Signac)
library(duckplyr)
library(cowplot)
library(TRIPOD)

prep_path = here(
    "processed-data", "13_tripod_trios", "02_tripod_preprocess",
    "preprocessed_objects_%s.qs2"
)
trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
model_paths = here(
    "processed-data", "13_tripod_trios", "03_tripod_trios", "fit_models_%s.qs2"
)
plot_dir = here("plots", "13_tripod_trios", "10_scatter_checks")
cell_types = c("MHb.2", "LHb.2.7")
num_examples = 5

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

trio_df = read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    filter(cell_type %in% cell_types, is_intersect, is_unique, is_top_TF) |>
    collect() |>
    group_by(cell_type, peak, gene, TF) |>
    filter(n() == 2) |>
    ungroup() |>
    distinct(cell_type, peak, gene, TF)

#-------------------------------------------------------------------------------
#   Scatter plots: expression vs. accessibility, colored by TF expression
#-------------------------------------------------------------------------------

xy_mat_list = lapply(
    unique(trio_df$cell_type),
    function(this_cell_type) {
        model_paths |>
            sprintf(this_cell_type) |>
            qs_read()
    }
)
names(xy_mat_list) = unique(trio_df$cell_type)

plot_list = list()
for (ct in unique(trio_df$cell_type)) {
    mats = qs_read(sprintf(prep_path, ct))$metacell_seur
    ct_trios = trio_df |> filter(cell_type == ct)

    for (i in seq_len(nrow(ct_trios))) {
        row    = ct_trios[i, ]
        expr   = mats$rna[, row$gene, drop = TRUE]
        access = mats$peak[, row$peak, drop = TRUE]
        tf_expr = mats$rna[, row$TF, drop = TRUE]

        plot_list[[length(plot_list) + 1]] = tibble(
                cell_type   = ct,
                gene        = row$gene,
                TF          = row$TF,
                access      = access,
                expr        = expr,
                tf_expr_log = log2(tf_expr + 1)
            ) |>
            ggplot(aes(x = access, y = expr, color = tf_expr_log)) +
                geom_point() +
                scale_color_viridis_c() +
                labs(
                    title = sprintf(
                        "%s\n%s | %s", ct, trio_df$gene[i], trio_df$TF[i]
                    ),
                    x = "Access.", y = "Expr."
                ) +
                theme_bw(base_size = 15) +
                theme(legend.position = "none", plot.title = element_text(size = 10))
      
        # for (this_test_level in 1:2) {
        #     plot_list[[length(plot_list) + 1]] = plotGenePeakTFScatter(
        #         xymats = xy_mat_list[[ct]][[ct_trios$gene[i]]],
        #         peak.name = ct_trios$peak[i], TF.name = ct_trios$TF[i],
        #         to.plot = "TRIPOD", match.by = "Yj",
        #         level = this_test_level, cap.at.quantile = 0
        #     )
        # }
    }
}

pdf(file.path(plot_dir, "top_trio_scatter.pdf"))
plot_grid(plotlist = plot_list, ncol = 3)
dev.off()

session_info()
  