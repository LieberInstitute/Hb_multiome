library(tidyverse)
library(here)
library(ComplexHeatmap)
library(RColorBrewer)
library(sessioninfo)
library(getopt)

# Import command-line parameters
spec <- matrix(
    c(
        c("gwas", "cell_type_group"),
        c("g", "c"),
        rep("1", 2),
        rep("character", 2),
        rep("Add variable description here", 2)
    ),
    ncol = 5
)
opt <- getopt(spec)

message("Using the following parameters:")
print(opt)

results_path = here(
    'processed-data', '10_MAGMA', opt$gwas,
    sprintf('%s.gsa.out', opt$cell_type_group)
)
gene_sets_path = here(
    'processed-data', '10_MAGMA', 'gene_sets',
    sprintf('%s.tsv', opt$cell_type_group)
)
gene_stat_path = here(
    'processed-data', '10_MAGMA', opt$gwas,
    sprintf('%s.genes.out', opt$gwas)
)
plot_dir = here('plots', '10_MAGMA', opt$gwas)
sig_cutoff = 0.05

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

#   MAGMA outputs have a variable amount of header lines. Auto-detect then
#   header length and read in dynamically
read_table_auto_skip = function(path, check_lines = 100) {
    n_skip = sum(grepl('^#', readLines(path, n = check_lines)))
    clean_df = read_table(path, skip = n_skip, show_col_types = FALSE)
    return(clean_df)
}

#   Read in and clean MAGMA results
results_df = read_table_auto_skip(results_path) |>
    mutate(
        cell_type = str_extract(FULL_NAME, '^[^_]+'),
        peak_category = str_extract(FULL_NAME, '(?<=_).+'),
        neg_log_p = -log10(P)
    ) |>
    select(cell_type, peak_category, neg_log_p)

#   Read in the actual gene sets, since we'll later count the size of the union
#   of genes across cell types and peak categories. Here we only count genes
#   with statistics from MAGMA
gene_stat_df = read_table(gene_stat_path, show_col_types = FALSE)
gene_sets_df = read_tsv(gene_sets_path, show_col_types = FALSE) |>
    mutate(
        cell_type = str_extract(set_id, '^[^_]+'),
        peak_category = str_extract(set_id, '(?<=_).+')
    ) |>
    select(-set_id) |>
    filter(
        cell_type %in% results_df$cell_type,
        peak_category %in% results_df$peak_category,
        gene_id %in% gene_stat_df$GENE
    )

#   Number of genes per cell type (note sorting ensure proper ordering of cell
#   types)
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

#   Number of genes per peak category (note sorting ensure proper ordering of
#   peak categories)
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

#   Form a matrix of p-values, filling in NAs for missing combos of cell type
#   and peak category
results_mat = results_df |>
    arrange(cell_type, peak_category) |>
    pivot_wider(
        names_from = peak_category,
        values_from = neg_log_p,
        values_fill = NA_real_
    ) |>
    column_to_rownames('cell_type') |>
    as.matrix()

#   Ensure ordering matches row and column annotations
results_mat = results_mat[
    sort(rownames(results_mat)), sort(colnames(results_mat))
]

#   Annotate significant p-values only
anno_mat = ifelse(
    !is.na(results_mat) & (results_mat >= -log10(sig_cutoff)),
    round(results_mat, 2),
    ''
)

pdf(
    file.path(plot_dir, sprintf('%s_pval_heatmap.pdf', opt$cell_type_group)),
    width = 5
)
Heatmap(
    results_mat,
    name = '-log10(p)',
    row_title = 'Cell Type',
    column_title = 'Peak Category',
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    top_annotation = col_anno,
    right_annotation = row_anno,
    col = c("white", colorRampPalette(brewer.pal(9, "YlOrRd"))(50)),
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(anno_mat[i, j], x, y, gp = gpar(fontsize = 10))
    }
)
dev.off()

session_info()
