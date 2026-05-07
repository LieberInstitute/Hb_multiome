library(tidyverse)
library(Seurat)
library(Signac)
library(here)
library(qs2)
library(sessioninfo)

seur_path = here(
    "processed-data", "11_link_prep", "05_metacell_aggregate",
    "meta_seur.qs2"
)
plot_dir = here("plots", "11_link_prep", "06_metacell_summaries")
cell_type_colors = c(
    "Oligo" = "#4d5802",
    "OPC" = "#d3c871",
    "Microglia" = "#222222",
    "Astrocyte" = "#8d363c",
    "Endo" = "#ee6c14",
    "Inhib.Thal" = '#b5a2ff',
    "Excit.Thal" = "#9e4ad1",
    "LHb.1.3.4" = "#0085af",
    "LHb.2.7" = "#7DF9FF",
    "LHb.4" = '#76c2af',
    "MHb.1" = "#FF00FF",
    "MHb.1.2" = "#c76a6a",
    "MHb.2" = "#FAA0A0",
    "MHb.3" = "#fa246a",
    "Inhib_LHb_4.1" = "#23008b",
    "Inhib_LHb_4.2" = "#3333fa",
    "Ependymal" = "#f5a105ff"
)

dir.create(plot_dir, showWarnings = FALSE)

seur = qs_read(seur_path)

meta_df = seur@meta.data |>
    as_tibble() |>
    select(refined_mid_cluster, size)

## Order cell types by median size for violin plot
violin_order = meta_df |>
    summarize(med = median(size), .by = refined_mid_cluster) |>
    arrange(med) |>
    pull(refined_mid_cluster)

## Order cell types by descending count for barplot
bar_order = meta_df |>
    count(refined_mid_cluster) |>
    arrange(desc(n)) |>
    pull(refined_mid_cluster)

## Violin plot: metacell size by cell type
pdf(file.path(plot_dir, "violin_metacell_size_by_celltype.pdf"), width = 10, height = 6)
meta_df |>
    mutate(refined_mid_cluster = factor(refined_mid_cluster, levels = violin_order)) |>
    ggplot(aes(x = refined_mid_cluster, y = size, fill = refined_mid_cluster)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.7) +
    scale_fill_manual(values = cell_type_colors) +
    labs(
        title = "Metacell Size by Cell Type",
        x = "Cell Type",
        y = "Metacell Size (# cells)"
    ) +
    theme_bw(base_size = 15) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
dev.off()

## Bar plot: number of metacells per cell type
pdf(file.path(plot_dir, "barplot_metacells_per_celltype.pdf"), width = 10, height = 6)
meta_df |>
    count(refined_mid_cluster) |>
    mutate(refined_mid_cluster = factor(refined_mid_cluster, levels = bar_order)) |>
    ggplot(aes(x = refined_mid_cluster, y = n, fill = refined_mid_cluster)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.4, size = 4) +
    scale_fill_manual(values = cell_type_colors) +
    labs(
        title = "Number of Metacells per Cell Type",
        x = "Cell Type",
        y = "Number of Metacells"
    ) +
    theme_bw(base_size = 15) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
    )
dev.off()

message('MHb.3, despite having as many metacells as there are donors, has a mix of donors per metacell:')
seur@meta.data |>
    as_tibble() |>
    select(refined_mid_cluster, orig.ident_purity) |>
    filter(refined_mid_cluster == 'MHb.3') |>
    print(n = 10)

session_info()
