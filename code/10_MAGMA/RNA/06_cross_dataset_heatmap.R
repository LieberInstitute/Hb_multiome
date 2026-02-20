library(tidyverse)
library(here)
library(viridis)
library(sessioninfo)

multiome_path = here('processed-data', '10_MAGMA', 'RNA', 'heatmap_results.csv')
extra_hd_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/no_secondary/MAGMA/extracellular/heatmap_data.csv'
cell_hd_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/no_secondary/MAGMA/heatmap_data.csv'
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

################################################################################
#   Functions
################################################################################

p_val_heatmap = function(results_df, gwas_groups, f_name) {
    p = ggplot(
            results_df,
            aes(
                x = cell_type, y = dataset, fill = neg_log_p, label = p_label
            )
        ) +
        geom_tile() +
        geom_text(size = 6) +
        scale_fill_viridis_c() +
        facet_wrap(~gwas_group, nrow = 1) +
        theme_bw(base_size = 20) +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text.y.right = element_text(angle = 0)
        ) +
        labs(x = "Cell Type", y = "GWAS Trait", fill = "-log10(p)")
    pdf(
        file.path(plot_dir, f_name),
        width = 3 + 0.3 * length(gwas_groups) * length(unique(results_df$cell_type)),
        height = 5
    )
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

multiome_df = read_csv(multiome_path, show_col_types = FALSE) |>
    mutate(dataset = "multiome")
extra_hd_df = read_csv(extra_hd_path, show_col_types = FALSE) |>
    mutate(dataset = "extra_HD")
cell_hd_df = read_csv(cell_hd_path, show_col_types = FALSE) |>
    mutate(dataset = "cell_HD")
results_df = bind_rows(multiome_df, extra_hd_df, cell_hd_df)

#   P-value heatmaps split by substance-use-related traits vs. others
for (gwas_set in names(gwas_groups)) {
    for (cell_type_group in cell_type_groups) {
        p_val_heatmap(
            results_df = results_df |>
                filter(
                    gwas_group %in% gwas_renaming[gwas_groups[[gwas_set]]],
                    cell_type_group == !!cell_type_group
                ),
            gwas_groups = gwas_renaming[gwas_groups[[gwas_set]]],
            f_name = sprintf(
                "cross_dataset_%s_%s.pdf", gwas_set, cell_type_group
            )
        )
    }
}

#   P-value heatmap for the 5 factors + P factor
for (cell_type_group in cell_type_groups) {
    p_val_heatmap(
        results_df = results_df |>
            filter(!is.na(gwas_factor), cell_type_group == !!cell_type_group) |>
            mutate(gwas_group = gwas_factor),
        gwas_groups = gwas_factors,
        f_name = sprintf("cross_dataset_5_factors_%s.pdf", cell_type_group)
    )
}

session_info()
