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
gwas_groups = list(
    substance = c(
        'SUD2020', 'AUD', 'CUD', 'ext_cannabis', 'lifetime_cannabis', 'OUD',
        'SUD2', 'SUD3'
    ),
    non_substance = c(
        'MDD2019', 'panic', 'SCZ', 'compulsive', 'internalizing', 'neurodev',
        'p_factor', 'SCZ_BPD'
    )
)
gwas_renaming = c(
    'MDD2019' = 'MDD',
    'panic' = 'Panic Disorder',
    'compulsive' = 'Compuls. Dis.',
    'SCZ' = 'SCZ',
    'SCZ_BPD' = 'SCZ/BPD',
    'AUD' = 'AUD',
    'CUD' = 'CUD',
    'ext_cannabis' = 'Ext. Cannabis',
    'lifetime_cannabis' = 'Life. Cannabis',
    'SUD2020' = 'OUD 1',
    'OUD' = 'OUD 2',
    'SUD2' = 'SUD 1',
    'SUD3' = 'SUD 2',
    'internalizing' = 'Intern. Disorders',
    'neurodev' = 'Neurodev.',
    'p_factor' = 'P Factor'
)
peak_levels = c(
    'Linked_DAR_enriched', 'Linked_DAR_depleted', 'Linked_DAR_discordant',
    'Linked_OCR_enriched', 'Linked_OCR_depleted'
)
#   This isn't respected in facet_grid(scales = 'free_y')
# cell_type_levels = c(
#     'Hb', 'MHb', 'LHb', 'MHb.1', 'MHb.2', 'MHb.1.2', 'MHb.3', 'LHb.1',
#     'LHb.1.3', 'LHb.1.3.4', 'LHb.4', 'LHb.2.7', 'Astrocyte', 'Endo',
#     'Excit.Thal', 'Inhib.Thal', 'Thal', 'Microglia', 'Oligo', 'OPC'
# )
sig_cutoff = 0.05

################################################################################
#   Functions
################################################################################

#   Heatmap grid where rows are cell-type groups and columns are GWAS groups
p_val_heatmap = function(results_df, gwas_groups, f_name) {
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
                cell_type_group ~ gwas_group, scales = "free_y",
                space = "free_y"
            ) +
            theme_bw(base_size = 20) +
            theme(
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                strip.text.y.right = element_text(angle = 0)
            ) +
            labs(x = "Peak Category", y = "Cell Type", fill = "-log10(P)")
    pdf(
        file.path(plot_dir, f_name),
        width = 2 + 2 * length(gwas_groups),
        height = 3 + 2 * length(cell_type_groups)
    )
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

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
    for (gwas_group in unlist(gwas_groups)) {
        #   Read in and clean MAGMA results
        results_df_list[[length(results_df_list) + 1]] = read_table_auto_skip(
                sprintf(results_path, gwas_group, cell_type_group)
            ) |>
            mutate(
                cell_type = str_extract(FULL_NAME, '^[^_]+'),
                peak_category = factor(
                    str_extract(FULL_NAME, '(?<=_).+'), levels = peak_levels
                ),
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
        gwas_group = factor(gwas_renaming[gwas_group], levels = gwas_renaming),
        cell_type_group = factor(cell_type_group, levels = cell_type_groups)
    ) |>
    #   Remove redundant rows
    filter(grepl('Hb', cell_type) | cell_type_group == 'mid')

stopifnot(!any(is.na(results_df$peak_category)))

#   P-value heatmaps split by substance-use-related traits vs. others
for (gwas_set in names(gwas_groups)) {
    p_val_heatmap(
        results_df = results_df |>
            filter(gwas_group %in% gwas_renaming[gwas_groups[[gwas_set]]]),
        gwas_groups = gwas_renaming[gwas_groups[[gwas_set]]],
        f_name = sprintf("shared_heatmap_%s.pdf", gwas_set)
    )
}

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

#   Fraction of significant peak categories per cell type and GWAS
p = results_df |>
    group_by(cell_type, cell_type_group, gwas_group) |>
    summarize(
        frac_signif = sum(neg_log_p > -log10(sig_cutoff)) /
            length(unique(results_df$peak_category))
    ) |>
    ggplot(aes(x = gwas_group, y = cell_type, fill = frac_signif)) +
        geom_tile() +
        scale_fill_viridis_c() +
        facet_grid(
            rows = vars(cell_type_group), scales = "free_y", space  = "free_y"
        ) +
        theme_bw(base_size = 20) +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text.y.right = element_text(angle = 0)
        ) +
        labs(x = "GWAS", y = "Cell Type", fill = "Fraction\nSignificant")
pdf(
    file = file.path(plot_dir, "shared_heatmap_substance_frac_signif.pdf"),
    width = 9, height = 6
)
print(p)
dev.off()

session_info()
  