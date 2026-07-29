library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(qs2)
library(ComplexHeatmap)
library(circlize)

cell_map_path = here("raw-data", "cell_type_map.csv")
colors_path = here('code', '05_03_annotation_adjustments', 'celltype_colors.R')
dar_path = here(
    'processed-data', '15_DARs', '07_cell_level_gather', 'DARs_all.csv.gz'
)
seur_path = here(
    'processed-data', '11_link_prep', '02_rebuild_atac_assay',
    'cell_level_seur.qs2'
)
out_path = here('plots', '15_DARs', '05_heatmap')
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

source(colors_path)
cell_type_colors = my_colors_mid

cluster_map = read_csv(cell_map_path, show_col_types = FALSE)
rename_map = stats::setNames(
    cluster_map$new_cell_type, cluster_map$old_cell_type
)

seur = qs_read(seur_path)
dar_df = read_csv(dar_path, show_col_types = FALSE) |>
    filter(avg_log2FC > 0)

## Top 5 peaks per cell type by FC ----------------------------

top_peaks = dar_df |>
    group_by(cell_type) |>
    slice_max(avg_log2FC, n = 5, with_ties = FALSE) |>
    ungroup()

peak_vec = unique(top_peaks$peak)
cell_types_ordered = unique(top_peaks$cell_type)

# Translate old cell type names to new display names, ordered as in cluster_map
cell_types_display = dplyr::coalesce(
    unname(rename_map[cell_types_ordered]),
    cell_types_ordered
)
cell_types_display = factor(cell_types_display, levels = cluster_map$new_cell_type)
cell_types_display = levels(droplevels(cell_types_display))

## Extract ATAC data layer, average by cell type, Z-score rows --------------

mat = LayerData(seur, assay = "ATAC", layer = "data")[peak_vec, ]
meta = seur@meta.data

cell_type_means = sapply(cell_types_ordered, function(ct) {
    cols = rownames(meta)[meta$refined_mid_cluster == ct]
    if (length(cols) == 1) mat[, cols]
    else rowMeans(mat[, cols])
})

# Z-score each peak (row) across cell types
mat_z = t(scale(t(cell_type_means)))

## Order peaks by canonical cell type order, then by p_val_adj within group --

peak_to_ct = top_peaks |>
    select(peak, cell_type) |>
    distinct() |>
    group_by(peak) |>
    slice(1) |>
    ungroup() |>
    mutate(
        cell_type_display = dplyr::coalesce(
            unname(rename_map[cell_type]), cell_type
        ) |>
            factor(levels = cell_types_display)
    )

peak_order = top_peaks |>
    mutate(
        cell_type_display = dplyr::coalesce(
            unname(rename_map[cell_type]), cell_type
        ) |>
            factor(levels = cell_types_display)
    ) |>
    arrange(cell_type_display, p_val_adj) |>
    pull(peak) |>
    unique()

row_ct = peak_to_ct$cell_type_display[match(peak_order, peak_to_ct$peak)]

## Build heatmap (x = peaks, y = cell types) --------------------------------

ct_colors = cell_type_colors[cell_types_display]

col_fun = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

# Transpose: rows = cell types, columns = peaks
# Row order matches the column grouping order so the diagonal is aligned
mat_t = t(mat_z[peak_order, ])
# mat_z columns are old cell type names; rename rows to display names for the heatmap
rownames(mat_t) = dplyr::coalesce(unname(rename_map[rownames(mat_t)]), rownames(mat_t))
mat_t = mat_t[cell_types_display, ]

col_split = factor(row_ct, levels = cell_types_display)

col_ha = columnAnnotation(
    cell_type = row_ct,
    col = list(cell_type = ct_colors),
    annotation_name_gp = gpar(fontsize = 9),
    show_legend = FALSE
)

ht = Heatmap(
    mat_t,
    name = "Z-score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 90,
    row_names_gp = gpar(fontsize = 9),
    row_names_side = "left",
    column_title = "Top DAR Peaks per Cell Type",
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    top_annotation = col_ha,
    column_split = col_split,
    column_gap = unit(1, "mm"),
    border = TRUE,
    heatmap_legend_param = list(title = "Z-score", title_gp = gpar(fontsize = 9))
)

## Save ---------------------------------------------------------------------

pdf(file.path(out_path, "dar_accessibility_heatmap.pdf"), width = 9, height = 6)
draw(ht)
dev.off()

## Session info -------------------------------------------------------------

session_info()
