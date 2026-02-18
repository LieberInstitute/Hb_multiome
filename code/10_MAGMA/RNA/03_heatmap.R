library(tidyverse)
library(here)
library(viridis)
library(sessioninfo)

results_path = here('processed-data','10_MAGMA', 'RNA', '%s', '%s.gsa.out')
plot_dir = here('plots', '10_MAGMA', 'RNA')

cell_type_groups = c('broad', 'fine')
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

#   Heatmap grid where rows are cell-type groups and columns are GWAS groups
p_val_heatmap = function(results_df, gwas_groups, f_name) {
    df_plot = results_df |>
    filter(gwas_group %in% gwas_groups) |>
    mutate(gwas_group = factor(gwas_group, levels = gwas_groups))

  p = df_plot |>
    ggplot(aes(x = gwas_group, y = cell_type, fill = neg_log_p, label = p_label)) +
      geom_tile() +
      geom_text(size = 6) +
      scale_fill_viridis_c() +
      facet_grid(
        rows = vars(cell_type_group),
        scales = "free_y",
        space  = "free_y"
      ) +
      theme_bw(base_size = 20) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.y.right = element_text(angle = 0)
      ) +
      labs(x = "GWAS", y = "Cell Type", fill = "-log10(P)")

  pdf(
    file.path(plot_dir, f_name),
    width  = 2 + 0.9 * length(gwas_groups),
    height = 3 + 3 * length(cell_type_groups)
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

    results_df_list[[length(results_df_list) + 1]] =
      read_table_auto_skip(sprintf(results_path, gwas_group, cell_type_group)) |>
      mutate(
        cell_type = VARIABLE,
        neg_log_p = -log10(pmax(P, .Machine$double.xmin)),
        cell_type_group = cell_type_group,
        gwas_group = gwas_group
      ) |>
      select(
        cell_type, neg_log_p, cell_type_group, gwas_group
      )
  }
}

#   Gather into one tibble
results_df = bind_rows(results_df_list) |>
  mutate(
    p_label = ifelse(neg_log_p > -log10(sig_cutoff), "*", ""),
    is_signif = neg_log_p > -log10(sig_cutoff),
    gwas_group = factor(gwas_renaming[gwas_group], levels = gwas_renaming),
    cell_type_group = factor(cell_type_group, levels = cell_type_groups)
  )
  
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
  select(cell_type, neg_log_p, gwas_group, cell_type_group) |>
  pivot_wider(names_from = gwas_group, values_from = neg_log_p) |>
  filter(!is.na(`OUD 1`) & !is.na(`OUD 2`)) |>
  group_by(cell_type_group) |>
  summarise(OUD_cor = cor(`OUD 1`, `OUD 2`)) |>
  print()

#   Fraction of significant peak categories per cell type and GWAS
p = results_df |>
  group_by(cell_type, cell_type_group, gwas_group) |>
  summarize(is_signif = any(is_signif), .groups = "drop") |>
  ggplot(aes(x = gwas_group, y = cell_type, fill = as.numeric(is_signif))) +
    geom_tile() +
    scale_fill_viridis_c() +
    facet_grid(
      rows = vars(cell_type_group), scales = "free_y", space = "free_y"
    ) +
    theme_bw(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text.y.right = element_text(angle = 0)
    ) +
    labs(x = "GWAS", y = "Cell Type", fill = "Significant\n(0/1)")

pdf(
  file = file.path(plot_dir, "shared_heatmap_substance_signif.pdf"),
  width = 9, height = 8
)
print(p)
dev.off()

session_info()