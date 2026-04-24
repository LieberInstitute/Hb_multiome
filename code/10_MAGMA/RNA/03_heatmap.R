library(tidyverse)
library(here)
library(viridis)
library(sessioninfo)

results_path = here(
    'processed-data','10_MAGMA', 'RNA', '%s', 'mean_ratio','%s.gsa.out'
)
gwas_name_path = here('processed-data', '10_MAGMA', 'RNA', 'gwas_info.csv')
low_genes_path = here('processed-data', '10_MAGMA', 'RNA', 'low_gene_sets.csv')
out_path = here('processed-data', '10_MAGMA', 'RNA', 'heatmap_data.csv')
plot_dir = here('plots', '10_MAGMA', 'RNA')

cell_type_groups = c('broad', 'mid', 'fine')
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

gwas_map = read_csv(gwas_name_path, show_col_types = FALSE) |>
    mutate(nickname = ifelse(nickname == 'MDD2019', 'MDD', nickname)) |>
    select(nickname, manuscript_name)

#   Read in all GWAS results for all cell-type groups
results_df_list = list()
for (cell_type_group in cell_type_groups) {
    for (gwas_group in gwas_map$nickname) {
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

low_genes_df = read_csv(low_genes_path, show_col_types = FALSE) |>
    dplyr::rename(cell_type_group = cell_type_res) |>
    mutate(low_gene_count = TRUE)

results_df = bind_rows(results_df_list) |>
    left_join(gwas_map, by = c('gwas_group' = 'nickname')) |>
    left_join(
        low_genes_df, by = c('cell_type_group', 'cell_type', 'manuscript_name')
    ) |>
    mutate(
        cell_type_group = factor(cell_type_group, levels = cell_type_groups),
        low_gene_count = ifelse(is.na(low_gene_count), FALSE, low_gene_count),
        p_label = ifelse(neg_log_p > -log10(sig_cutoff), "*", ""),
        #   Add a question mark for results determined from small gene sets
        p_label = ifelse(low_gene_count, paste0(p_label, "?"), p_label)
    ) |>
    dplyr::rename(gwas_nickname = gwas_group, gwas_group = manuscript_name)

write_csv(results_df, out_path)

#   P-value heatmaps split by substance-use-related traits vs. others
gwas_categories = list(
    substance = gwas_map$manuscript_name[
        grepl('^[ACSO]UD_', gwas_map$manuscript_name)
    ],
    psychiatric = gwas_map$manuscript_name[
        !grepl('^[ACSO]UD_|^p_factor_', gwas_map$manuscript_name)
    ],
    factor = c(
        gwas_map$manuscript_name[grepl('_F[1-5]_', gwas_map$manuscript_name)] |>
            (\(x) x[order(as.integer(stringr::str_extract(x, '(?<=_F)\\d(?=_)')))])(),
        'p_factor_Grotzinger'
    )
)

for (gwas_category in names(gwas_categories)) {
    p_val_heatmap(
        results_df = results_df |>
            filter(gwas_group %in% gwas_categories[[gwas_category]]) |>
            mutate(
                gwas_group = factor(
                    gwas_group, levels = gwas_categories[[gwas_category]]
                )
            ),
        gwas_groups = gwas_categories[[gwas_category]],
        f_name = sprintf("heatmap_%s.pdf", gwas_category)
    )
}

session_info()
