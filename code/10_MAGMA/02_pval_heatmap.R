library(tidyverse)
library(here)
library(ComplexHeatmap)
library(RColorBrewer)
library(sessioninfo)

results_path = here(
    'processed-data', '10_MAGMA', 'first_test', 'first_test.gsa.out'
)
gene_sets_path = here(
    'processed-data', '10_MAGMA', 'first_test', 'input_gene_sets.tsv'
)
plot_dir = here('plots', '10_MAGMA')

dir.create(plot_dir, showWarnings = FALSE)

results_df = read_table(results_path, skip = 3, show_col_types = FALSE) |>
    mutate(
        cell_type = str_extract(FULL_NAME, '^[^_]+'),
        peak_category = str_extract(FULL_NAME, '(?<=_).+'),
        neg_log_p = -log10(P)
    ) |>
    select(cell_type, peak_category, neg_log_p)

gene_sets_df = read_tsv(gene_sets_path, show_col_types = FALSE) |>
    mutate(
        cell_type = str_extract(set_id, '^[^_]+'),
        peak_category = str_extract(set_id, '(?<=_).+')
    ) |>
    select(-set_id) |>
    filter(
        cell_type %in% results_df$cell_type,
        peak_category %in% results_df$peak_category
    )

row_anno = HeatmapAnnotation(
    n_genes = anno_barplot(
        gene_sets_df |>
            group_by(cell_type) |>
            summarise(n = length(unique(gene_id))) |>
            arrange(cell_type) |>
            pull(n)
    ),
    which = 'row'
)

col_anno = HeatmapAnnotation(
    n_genes = anno_barplot(
        gene_sets_df |>
            group_by(peak_category) |>
            summarise(n = length(unique(gene_id))) |>
            arrange(peak_category) |>
            pull(n)
    ),
    which = 'column'
)

pdf(file.path(plot_dir, 'test_pval_heatmap.pdf'))
results_df |>
    arrange(cell_type, peak_category) |>
    pivot_wider(
        names_from = peak_category,
        values_from = neg_log_p,
        values_fill = NA_real_
    ) |>
    column_to_rownames('cell_type') |>
    as.matrix() |>
    Heatmap(
        name = '-log10(p)',
        row_title = 'Cell Type',
        column_title = 'Peak Category',
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        top_annotation = col_anno,
        right_annotation = row_anno,
        col = c("white", colorRampPalette(brewer.pal(9, "YlOrRd"))(50))
    )
dev.off()

session_info()
