library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(qs2)
library(sessioninfo)
library(ComplexHeatmap)
library(circlize)

seur_path = here(
    "processed-data", "05_Clustering_ARCr", "24_smaller_seurat",
    "general_purpose_seur.qs2"
)
marker_path = here(
    "processed-data", "10_MAGMA", "RNA", "gene_sets", "fine.tsv"
)
plot_dir = here("plots", "05_Clustering_ARCr", "25_manuscript_plots")
markers_per_cell_type = 4

dir.create(plot_dir, showWarnings = FALSE)

seur = qs_read(seur_path)

################################################################################
#   WNN UMAP at mid resolution
################################################################################

#   First, just a UMAP of the WNN results at mid resolution
p = DimPlot(
        seur,
        reduction = "wnn.umap",
        group.by = "mid_cluster",
        label = TRUE,
        label.size = 7,
        repel = TRUE
    ) +
    theme_bw(base_size = 20) &
    NoLegend()
png(file.path(plot_dir, "wnn_umap_mid.png"), width = 800, height = 800)
print(p)
dev.off()

################################################################################
#   Habenula-focused marker heatmap
################################################################################

#   Grab mid-resolution Hb markers (mean ratio)
marker_df = read_tsv(marker_path, show_col_types = FALSE) |>
    filter(grepl('^[ML]Hb', set_id)) |>
    group_by(set_id) |>
    slice_max(order_by = mean_ratio, n = markers_per_cell_type) |>
    ungroup() |>
    select(set_id, gene_name)

stopifnot(all(marker_df$gene_name %in% rownames(seur[['RNA']])))

#   Determine the average expression of these markers in each cell type
marker_df_list = list()
for (i in seq_len(nrow(marker_df))) {
    marker_df_list[[i]] = tibble(
        gene_name = marker_df$gene_name[i],
        marker_cell_type = marker_df$set_id[i],
        measured_cell_type = seur@meta.data$mid_cluster,
        value = unname(seur[['RNA']]$data[marker_df$gene_name[i], ])
    )
}
marker_df = bind_rows(marker_df_list) |>
    group_by(gene_name, marker_cell_type, measured_cell_type) |>
    summarize(value = mean(value)) |>
    ungroup()

# Prepare the matrix for the heatmap
heatmap_mat = marker_df |>
    select(gene_name, measured_cell_type, value) |>
    pivot_wider(names_from = measured_cell_type, values_from = value) |>
    column_to_rownames("gene_name") |>
    as.matrix()

# Get the marker cell type for each gene
gene_annotations = marker_df |>
    select(gene_name, marker_cell_type) |>
    distinct() |>
    mutate(region = ifelse(startsWith(marker_cell_type, "M"), "MHb", "LHb"))

# Order genes by their marker cell type (alphabetical)
ordered_cell_types = sort(colnames(heatmap_mat))
gene_order = gene_annotations |>
    mutate(marker_order = match(marker_cell_type, ordered_cell_types)) |>
    arrange(marker_order) |>
    pull(gene_name)

# Reorder matrix and annotations
heatmap_mat = heatmap_mat[gene_order, ordered_cell_types]
gene_annotations = gene_annotations |>
    slice(match(gene_order, gene_name))

# Z-score normalize by row
heatmap_mat_scaled = t(scale(t(heatmap_mat)))

# Create annotations and split
row_ha = rowAnnotation(
    Region = gene_annotations$region,
    col = list(Region = c("MHb" = "#9F7E69", "LHb" = "#0D316B"))
)
row_split = factor(gene_annotations$region, levels = c("MHb", "LHb"))

# Create heatmap
ht = Heatmap(
    heatmap_mat_scaled,
    name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = row_split,
    right_annotation = row_ha,
    row_title = NULL,
    column_names_rot = 90,
    heatmap_legend_param = list(title = "Expression\nZ-score")
)

pdf(file.path(plot_dir, "marker_heatmap.pdf"))
draw(ht)
dev.off()
