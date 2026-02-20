library(tidyverse)
library(here)
library(viridis)
library(sessioninfo)

results_path = here('processed-data','10_MAGMA', 'RNA', '%s', '%s.gsa.out')
out_path = here('processed-data', '10_MAGMA', 'RNA', 'heatmap_results.csv')
plot_dir = here('plots', '10_MAGMA', 'RNA')

cell_type_groups = c('broad', 'mid', 'fine')
gwas_groups = list(
    substance = c(
        'SUD2020', 'AUD', 'CUD', 'ext_cannabis', 'lifetime_cannabis', 'OUD',
        'SUD2', 'SUD3'
    ),
    non_substance = c(
        'MDD', 'panic', 'SCZ', 'compulsive', 'internalizing', 'neurodev',
        'p_factor', 'SCZ_BPD'
    )
)
#   The 5 factors + general P factor from the paper:
#   https://doi.org/10.1038/s41586-025-09820-3
gwas_factors = c(
    'compulsive' = 'F1: compulsive',
    'SCZ_BPD' = 'F2: SCZ/BPD',
    'neurodev' = 'F3: neurodev',
    'internalizing' = 'F4: intern.',
    'SUD3' = 'F5: SUD',
    'p_factor' = 'P Factor'
)
gwas_renaming = c(
    'MDD' = 'MDD',
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
sig_cutoff = 0.05

################################################################################
#   Functions
################################################################################

p_val_heatmap = function(results_df, gwas_groups, f_name) {
    p = ggplot(
            results_df,
            aes(
                x = gwas_group, y = cell_type, fill = neg_log_p, label = p_label
            )
        ) +
        geom_tile() +
        geom_text(size = 6) +
        scale_fill_viridis_c() +
        facet_wrap(~cell_type_group, ncol = 3, scales = "free_y") +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(x = "GWAS Trait", y = "Cell Type", fill = "-log10(p)")
    pdf(
        file.path(plot_dir, f_name),
        width = 3 + 2 * length(gwas_groups),
        height = 3 + length(cell_type_groups)
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
            dplyr::rename(cell_type = VARIABLE) |>
            mutate(
                neg_log_p = -log10(P),
                cell_type_group = cell_type_group,
                gwas_group = gwas_group
            ) |>
            select(cell_type, neg_log_p, cell_type_group, gwas_group)
    }
}

#   Gather into one tibble
results_df = bind_rows(results_df_list) |>
    mutate(
        p_label = ifelse(neg_log_p > -log10(sig_cutoff), "*", ""),
        gwas_factor = factor(
            ifelse(
                gwas_group %in% names(gwas_factors),
                gwas_factors[gwas_group], NA
            ),
            levels = gwas_factors
        ),
        gwas_group = factor(gwas_renaming[gwas_group], levels = gwas_renaming),
        cell_type_group = factor(cell_type_group, levels = cell_type_groups)
    )

write_csv(results_df, out_path)

#   P-value heatmaps split by substance-use-related traits vs. others
for (gwas_set in names(gwas_groups)) {
    p_val_heatmap(
        results_df = results_df |>
            filter(gwas_group %in% gwas_renaming[gwas_groups[[gwas_set]]]),
        gwas_groups = gwas_renaming[gwas_groups[[gwas_set]]],
        f_name = sprintf("heatmap_%s.pdf", gwas_set)
    )
}

#   P-value heatmap for the 5 factors + P factor
p_val_heatmap(
    results_df = results_df |>
        filter(!is.na(gwas_factor)) |>
        mutate(gwas_group = gwas_factor),
    gwas_groups = gwas_factors,
    f_name = "heatmap_5_factors.pdf"
)

session_info()

session_info()