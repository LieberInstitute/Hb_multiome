library(sessioninfo)
library(tidyverse)
library(here)
library(ComplexHeatmap)
library(circlize)
library(duckplyr)

trio_paths = here(
    "processed-data", "13_tripod_trios", "03_tripod_trios", "trios_%s.parquet"
)
plot_dir = here("plots", "13_tripod_trios", "08_trio_heatmap")
cell_type1 = "MHb.2"
cell_type2 = "LHb.2.7"
FDR_thres = 0.05

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

trio_df_list = list()
for (cell_type in c(cell_type1, cell_type2)) {
    level_1 = read_parquet_duckdb(
            sprintf(trio_paths, cell_type), prudence = 'lavish'
        ) |>
        filter(condition_on == 'Yj', stringency_level == 1, adj < FDR_thres) |>
        dplyr::rename(adj_1 = adj, cor_val = coef) |>
        select(peak, gene, TF, adj_1, cor_val)
    level_2 = read_parquet_duckdb(
            sprintf(trio_paths, cell_type), prudence = 'lavish'
        ) |>
        filter(condition_on == 'Yj', stringency_level == 2, adj < FDR_thres) |>
        dplyr::rename(adj_2 = adj) |>
        select(peak, gene, TF, adj_2, coef)
    trio_df_list[[cell_type]] = inner_join(
            level_1, level_2, by = c("peak", "gene", "TF")
        ) |>
        mutate(cell_type = cell_type)
}
trio_df = bind_rows(trio_df_list) |>
    collect()

#-------------------------------------------------------------------------------
#   Side-by-side heatmap: cor_val and coef per cell type
#-------------------------------------------------------------------------------

# Shared color ranges across both cell types
cor_range  <- range(trio_df$cor_val)
coef_range <- range(trio_df$coef)

# Viridis for cor_val, plasma for coef — both perceptually uniform and positive
col_cor  <- colorRamp2(
    seq(cor_range[1],  cor_range[2],  length.out = 9),
    viridisLite::viridis(9)
)
col_coef <- colorRamp2(
    seq(coef_range[1], coef_range[2], length.out = 9),
    viridisLite::rocket(9)
)

make_heatmap <- function(df, ct) {
    sub <- df |>
        dplyr::filter(cell_type == ct) |>
        dplyr::arrange(desc(cor_val)) |>
        dplyr::mutate(label = paste(gene, TF, sep = "\n"))

    cor_mat  <- matrix(sub$cor_val, ncol = 1,
        dimnames = list(sub$label, "cor_val\n(level 1 coef)"))
    coef_mat <- matrix(sub$coef,    ncol = 1,
        dimnames = list(sub$label, "coef\n(level 2)"))

    ht_cor <- Heatmap(
        cor_mat,
        col = col_cor,
        name = "cor_val",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_title = ct,
        column_title_gp = gpar(fontsize = 11, fontface = "bold"),
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 9),
        show_heatmap_legend = (ct == cell_type1),
        cell_fun = function(j, i, x, y, w, h, fill)
            grid.text(sprintf("%.2f", cor_mat[i, j]), x, y,
                gp = gpar(fontsize = 7, col = "white"))
    )

    ht_coef <- Heatmap(
        coef_mat,
        col = col_coef,
        name = "coef",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        column_names_gp = gpar(fontsize = 9),
        show_heatmap_legend = (ct == cell_type1),
        cell_fun = function(j, i, x, y, w, h, fill)
            grid.text(sprintf("%.2f", coef_mat[i, j]), x, y,
                gp = gpar(fontsize = 7, col = "white"))
    )

    ht_cor + ht_coef
}

ht1 <- make_heatmap(trio_df, cell_type1)
ht2 <- make_heatmap(trio_df, cell_type2)

pdf(file.path(plot_dir, "trio_heatmap.pdf"), width = 6, height = 4)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1))
draw(ht1, newpage = FALSE)
popViewport()
pushViewport(viewport(layout.pos.col = 2))
draw(ht2, newpage = FALSE)
popViewport()
dev.off()

session_info()
