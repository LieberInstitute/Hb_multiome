library(sessioninfo)
library(tidyverse)
library(here)
library(qs2)
library(Seurat)
library(Signac)

source(here("code", "05_03_annotation_adjustments", "celltype_colors.R"))

prep_path = here(
    "processed-data", "13_tripod_trios", "02_tripod_preprocess",
    "preprocessed_objects_%s.qs2"
)
cell_map_path = here('raw-data', 'cell_type_map.csv')
cell_types = c(
    "Astrocyte", "Ependymal", "Excit.Thal", "Inhib_LHb_4.1",
    "Inhib_LHb_4.2", "Inhib.Thal", "LHb.1.3.4", "LHb.2.7", "LHb.4", "MHb.1",
    "MHb.1.2", "MHb.2", "Microglia", "Oligo", "OPC"
)
plot_dir = here("plots", "13_tripod_trios", "07_metacell_counts")

dir.create(plot_dir, showWarnings = FALSE)

cell_map_df = read_csv(cell_map_path, show_col_types = FALSE)

meta_df_list = list()
for (this_cell_type in cell_types) {
    pre_list = qs_read(sprintf(prep_path, this_cell_type))
    
    meta_df_list[[this_cell_type]] = pre_list$seur@meta.data |>
        as_tibble() |>
        group_by(seurat_clusters, orig.ident) |>
        summarize(donor_counts = n()) |>
        group_by(seurat_clusters) |>
        summarize(
            n_cells = n(),
            donor_purity = max(donor_counts) / sum(donor_counts)
        ) |>
        mutate(cell_type = this_cell_type) |>
        select(cell_type, n_cells, donor_purity)
}
meta_df = bind_rows(meta_df_list) |>
    left_join(cell_map_df, by = c("cell_type" = "old_cell_type")) |>
    select(new_cell_type, n_cells, donor_purity) |>
    dplyr::rename(cell_type = new_cell_type)

## Order cell types by median n_cells for violin plot
violin_order = meta_df |>
    summarize(med = median(n_cells), .by = cell_type) |>
    arrange(med) |>
    pull(cell_type)

## Order cell types by descending metacell count for barplot
bar_order = meta_df |>
    count(cell_type) |>
    arrange(desc(n)) |>
    pull(cell_type)

## Violin plot: metacell size by cell type
pdf(
    file.path(plot_dir, "violin_metacell_size_by_celltype.pdf"),
    width = 10, height = 4
)
meta_df |>
    mutate(cell_type = factor(cell_type, levels = violin_order)) |>
    ggplot(aes(x = cell_type, y = n_cells, fill = cell_type)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.7) +
    scale_fill_manual(values = my_colors_mid) +
    labs(x = "Cell Type", y = "Cells Per Metacell") +
    theme_bw(base_size = 16) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
dev.off()

## Bar plot: number of metacells per cell type
pdf(
    file.path(plot_dir, "barplot_metacells_per_celltype.pdf"),
    width = 10, height = 4
)
meta_df |>
    count(cell_type) |>
    mutate(cell_type = factor(cell_type, levels = bar_order)) |>
    ggplot(aes(x = cell_type, y = n, fill = cell_type)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.4, size = 6) +
    scale_fill_manual(values = my_colors_mid) +
    labs(x = "Cell Type", y = "Number of Metacells") +
    theme_bw(base_size = 16) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
dev.off()

## Order cell types by median donor_purity for violin plot
purity_order = meta_df |>
    summarize(med = median(donor_purity), .by = cell_type) |>
    arrange(med) |>
    pull(cell_type)

## Violin plot: donor purity by cell type
pdf(
    file.path(plot_dir, "violin_donor_purity_by_celltype.pdf"),
    width = 10, height = 4
)
meta_df |>
    mutate(cell_type = factor(cell_type, levels = purity_order)) |>
    ggplot(aes(x = cell_type, y = donor_purity, fill = cell_type)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.7) +
    scale_fill_manual(values = my_colors_mid) +
    labs(x = "Cell Type", y = "Donor Purity") +
    theme_bw(base_size = 16) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
dev.off()

session_info()
