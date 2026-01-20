#   There are too many combinations of results to show from 06_pval_heatmap. Do
#   a version where all GWAS sets and all cell-type resolutions are in a single
#   heatmap using ggplot2 and facet_grid rather than ComplexHeatmap

library(tidyverse)
library(here)
library(viridis)
library(sessioninfo)

results_path = here('processed-data', '10_MAGMA', '%s', '%s.gsa.out')
gene_sets_path = here('processed-data', '10_MAGMA', 'gene_sets', 'broad.tsv')
gene_stat_path = here('processed-data', '10_MAGMA', '%s', '%s.genes.out')
plot_dir = here('plots', '10_MAGMA')
cell_type_groups = c('broad', 'semi_broad', 'mid')
gwas_groups = c(
    'MDD', 'MDD2019', 'panic', 'SCZ', 'SUD2020', 'AUD', 'CUD', 'ext_cannabis',
    'lifetime_cannabis', 'OUD', 'SUD2', 'compulsive', 'internalizing',
    'neurodev', 'p_factor', 'SCZ_BPD', 'SUD3'
)
gwas_renaming = c(
    'MDD2019' = 'MDD',
    'panic' = 'Panic Disorder',
    'SCZ' = 'SCZ',
    'SUD2020' = 'OUD 1',
    'AUD' = 'AUD',
    'CUD' = 'CUD',
    'ext_cannabis' = 'Ext. Cannabis',
    'lifetime_cannabis' = 'Lifetime Cannabis',
    'OUD' = 'OUD 2',
    'SUD2' = 'SUD 1',
    'compulsive' = 'Compulsive Disorders',
    'internalizing' = 'Internalizing Disorders',
    'neurodev' = 'Neurodev. Disorders',
    'p_factor' = 'P Factor',
    'SCZ_BPD' = 'SCZ/BPD',
    'SUD3' = 'SUD 2'
)
sig_cutoff = 0.05

dir.create(plot_dir, showWarnings = FALSE)

#   MAGMA outputs have a variable amount of header lines. Auto-detect then
#   header length and read in dynamically
read_table_auto_skip = function(path, check_lines = 100) {
    n_skip = sum(grepl('^#', readLines(path, n = check_lines)))
    clean_df = read_table(path, skip = n_skip, show_col_types = FALSE)
    return(clean_df)
}

#   Read in all GWAS results for all cell-type groups
results_df_list = list()
for (cell_type_group in cell_type_groups) {
    for (gwas_group in gwas_groups) {
        #   Read in and clean MAGMA results
        results_df_list[[length(results_df_list) + 1]] = read_table_auto_skip(
                sprintf(results_path, gwas_group, cell_type_group)
            ) |>
            mutate(
                cell_type = str_extract(FULL_NAME, '^[^_]+'),
                peak_category = str_extract(FULL_NAME, '(?<=_).+'),
                neg_log_p = -log10(P),
                cell_type_group = cell_type_group,
                gwas_group = gwas_group
            ) |>
            select(
                cell_type, peak_category, neg_log_p, cell_type_group, gwas_group
            )
    }
}

#   Gather into one tibble
results_df = bind_rows(results_df_list) |>
    mutate(
        p_label = ifelse(neg_log_p > -log10(sig_cutoff), "*", ""),
        gwas_group = gwas_renaming[gwas_group],
        cell_type_group = factor(cell_type_group, levels = cell_type_groups)
    )

#   Heatmap grid where rows are cell-type groups and columns are GWAS groups
p = results_df |>
    ggplot(
            aes(
                x = peak_category, y = cell_type, fill = neg_log_p,
                label = p_label
            )
        ) +
        geom_tile() +
        geom_text(size = 6) +
        scale_fill_viridis_c() +
        facet_grid(
            cell_type_group ~ gwas_group, scales = "free_y", space = "free_y"
        ) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(x = "Peak Category", y = "Cell Type", fill = "-log10(P)")
pdf(
    file.path(plot_dir, "shared_heatmap.pdf"),
    width = 2 + 2 * length(gwas_groups), height = 3 + 2 * length(cell_type_groups)
)
print(p)
dev.off()

#   For each cell-type resolution, take the correlation of -log(p-value)
#   between the 2 OUD GWAS sets, which should be substantial and positive
message("Correlation of -log(p-value) between OUD 1 and OUD 2 for each cell-type resolution:")
results_df |>
    filter(gwas_group %in% c('OUD 1', 'OUD 2')) |>
    select(cell_type, peak_category, neg_log_p, gwas_group, cell_type_group) |>
    pivot_wider(names_from = gwas_group, values_from = neg_log_p) |>
    filter(!is.na(`OUD 1`) & !is.na(`OUD 2`)) |>
    group_by(cell_type_group) |>
    summarise(OUD_cor = cor(`OUD 1`, `OUD 2`)) |>
    print()

session_info()
