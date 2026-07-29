library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(qs2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

dar_path = here(
    'processed-data', '15_DARs', '07_cell_level_gather', 'DARs_all.csv.gz'
)
seur_path = here(
    'processed-data', '11_link_prep', '02_rebuild_atac_assay',
    'cell_level_seur.qs2'
)
out_path = here('plots', '15_DARs', '05_heatmap')
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

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

## Order rows by cell type grouping, then by p_val_adj within each group ----

peak_order = top_peaks |>
    arrange(cell_type, p_val_adj) |>
    pull(peak) |>
    unique()

mat_z_ordered = mat_z[peak_order, ]

peak_to_ct = top_peaks |>
    select(peak, cell_type) |>
    distinct() |>
    group_by(peak) |>
    slice(1) |>
    ungroup()

row_ct = peak_to_ct$cell_type[match(peak_order, peak_to_ct$peak)]

## Build heatmap (x = peaks, y = cell types) --------------------------------

n_ct = length(cell_types_ordered)
ct_colors = setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(n_ct),
    cell_types_ordered
)

col_fun = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

# Transpose: rows = cell types, columns = peaks
# Row order matches the column grouping order so the diagonal is aligned
mat_t = t(mat_z[peak_order, ])[cell_types_ordered, ]

col_split = factor(row_ct, levels = cell_types_ordered)

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
    column_names_gp = gpar(fontsize = 5),
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 8),
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

pdf(file.path(out_path, "dar_accessibility_heatmap.pdf"), width = 10, height = 14)
draw(ht)
dev.off()

## Session info -------------------------------------------------------------

session_info()
