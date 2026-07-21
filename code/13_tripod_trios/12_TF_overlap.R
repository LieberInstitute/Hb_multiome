library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(ComplexHeatmap)
library(circlize)

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
cell_map_path = here("raw-data", "cell_type_map.csv")
plot_dir = here("plots", "13_tripod_trios", "12_TF_overlap")

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

cluster_map = read_csv(cell_map_path, show_col_types = FALSE)
rename_map = stats::setNames(
    cluster_map$new_cell_type, cluster_map$old_cell_type
)

trio_df = read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    filter(is_intersect, stringency_level == 1) |>
    select(TF, cell_type) |>
    collect() |>
    mutate(
        cell_type = dplyr::coalesce(unname(rename_map[cell_type]), cell_type)
    )

#-------------------------------------------------------------------------------
#   TF breadth: histogram of how many cell types each TF appears in
#-------------------------------------------------------------------------------

p = trio_df |>
    distinct(TF, cell_type) |>
    count(TF, name = "n_cell_types") |>
    ggplot(aes(x = n_cell_types)) +
        geom_histogram(binwidth = 1) +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        scale_y_continuous(labels = scales::comma) +
        labs(x = "Number of cell types", y = "Number of TFs") +
        theme_bw(base_size = 15)
pdf(file.path(plot_dir, "TF_breadth_histogram.pdf"), width = 8, height = 6)
print(p)
dev.off()

#-------------------------------------------------------------------------------
#   TF set heatmaps
#-------------------------------------------------------------------------------

tf_sets = trio_df |>
    distinct(cell_type, TF) |>
    group_by(cell_type) |>
    summarise(tfs = list(TF), .groups = "drop")

ct_names = tf_sets$cell_type
n_ct = length(ct_names)

pair_stats = matrix(NA_real_, nrow = n_ct, ncol = n_ct,
                    dimnames = list(ct_names, ct_names))
jaccard_mat = pair_stats
overlap_mat = pair_stats

for (i in seq_len(n_ct)) {
    for (j in seq_len(n_ct)) {
        a = tf_sets$tfs[[i]]
        b = tf_sets$tfs[[j]]
        inter = length(intersect(a, b))
        jaccard_mat[i, j] = inter / length(union(a, b))
        overlap_mat[i, j] = inter / min(length(a), length(b))
    }
}

jaccard_range = range(jaccard_mat, na.rm = TRUE)
overlap_range = range(overlap_mat, na.rm = TRUE)

jaccard_col = colorRamp2(
    seq(jaccard_range[1], jaccard_range[2], length.out = 256),
    viridisLite::viridis(256)
)
overlap_col = colorRamp2(
    seq(overlap_range[1], overlap_range[2], length.out = 256),
    viridisLite::viridis(256)
)

# Jaccard index
ht_jaccard = Heatmap(
    jaccard_mat,
    name            = "Jaccard\nindex",
    col             = jaccard_col,
    cluster_rows    = TRUE,
    cluster_columns = TRUE,
    column_title    = "TF set Jaccard index between cell types",
    row_names_gp    = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 14),
    column_names_rot = 90
)

pdf(file.path(plot_dir, "Jaccard_TF_heatmap.pdf"), width = 9, height = 8)
draw(ht_jaccard)
dev.off()

# Intersection over smaller set
ht_overlap = Heatmap(
    overlap_mat,
    name            = "Overlap\ncoefficient",
    col             = overlap_col,
    cluster_rows    = TRUE,
    cluster_columns = TRUE,
    column_title    = "TF overlap coefficient between cell types",
    row_names_gp    = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 14),
    column_names_rot = 90
)

pdf(file.path(plot_dir, "Overlap_TF_heatmap.pdf"), width = 9, height = 8)
draw(ht_overlap)
dev.off()

#-------------------------------------------------------------------------------
#   Top 5 TFs per cell type by number of trios
#-------------------------------------------------------------------------------

p = trio_df |>
    count(cell_type, TF, name = "n_trios") |>
    group_by(cell_type) |>
    slice_max(n_trios, n = 5, with_ties = FALSE) |>
    ungroup() |>
    mutate(
        TF = paste(cell_type, TF, sep = "___"),
        TF = reorder(TF, n_trios)
    ) |>
    ggplot(aes(x = n_trios, y = TF)) +
    geom_col() +
    scale_y_discrete(labels = \(x) sub(".*___", "", x)) +
    facet_wrap(~cell_type, scales = "free") +
    labs(x = "Number of trios", y = NULL, title = "Top 5 TFs per cell type") +
    theme_bw(base_size = 13)

pdf(file.path(plot_dir, "Top_5_TFs.pdf"))
print(p)
dev.off()

session_info()
