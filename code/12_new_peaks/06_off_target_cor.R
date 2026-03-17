#   Using the data from every combination of two cell types, explore how the
#   distribution of correlation scores differ between the target and non-target
#   cell types. Also write a merged CSV of this data

library(here)
library(tidyverse)
library(sessioninfo)

cell_types = c(
    'Astrocyte', 'Endo', 'Excit.Thal', 'Inhib.Thal', 'LHb.1', 'LHb.1.3',
    'LHb.1.3.4', 'LHb.2.7', 'LHb.4', 'LHb.7', 'MHb.1', 'MHb.1.2', 'MHb.2',
    'MHb.3', 'Microglia', 'Oligo', 'OPC', 'Thal'
)
result_paths = here(
    'processed-data', '12_new_peaks', '01_link_peaks', '%s_%s.csv.gz'
)
out_path = here(
    'processed-data', '12_new_peaks', '01_link_peaks', 'all_data.csv.gz'
)
plot_dir = here("plots", "12_new_peaks")
fdr_cutoffs = c(0.1, 0.2, 0.3, 0.4, 0.425, 0.45, 0.475, 0.5)

#   First gather and write the key columns from every cell-type pair combination
result_df_list = list()
for (target_cell_type in cell_types) {
    for (other_cell_type in cell_types) {
        result_df_list[[other_cell_type]] = read_csv(
                sprintf(result_paths, target_cell_type, other_cell_type),
                show_col_types = FALSE
            ) |>
            select(peak, gene, score, FDR, target_cell_type, other_cell_type) |>
            mutate(
                target_cell_type = factor(
                    target_cell_type, levels = cell_types
                ),
                other_cell_type = factor(other_cell_type, levels = cell_types)
            )
    }
}
result_df = bind_rows(result_df_list)

write_csv(result_df, out_path)

result_df = result_df |>
    mutate(
        cell_type_group = ifelse(
            target_cell_type == other_cell_type, "target", "non_target"
        )
    ) |>
    select(cell_type_group, score, FDR)

#   Density plots of correlation scores split by whether target cell type matches
#   other cell type. No substantial difference is seen at this level
p = ggplot(
        result_df,
        aes(x = abs(score), fill = cell_type_group, color = cell_type_group)
    ) +
    geom_density(alpha = 0.5, linewidth = 1) +
    theme_bw(base_size = 20) +
    labs(
        x = "abs(Correlation Score)", y = "Density", fill = "Cell Type",
        color = "Cell Type"
    )
pdf(file.path(plot_dir, "off_target_cor_global.pdf"), width = 10, height = 5)
print(p)
dev.off()

result_df_list = list()
for (fdr_cutoff in fdr_cutoffs) {
    result_df_list[[as.character(fdr_cutoff)]] = result_df |>
        filter(FDR <= fdr_cutoff) |>
        mutate(FDR_class = sprintf("FDR <= %s", fdr_cutoff))
}
result_df = bind_rows(result_df_list) |>
    mutate(
        FDR_class = factor(
            FDR_class, levels = sprintf("FDR <= %s", fdr_cutoffs)
        )
    )

#   However when stratifying by FDR class, target-matching pairs have a
#   substantially larger magnitude of correlation for lower FDRs
p = ggplot(
        result_df,
        aes(x = abs(score), fill = cell_type_group, color = cell_type_group)
    ) +
    geom_density(alpha = 0.5, linewidth = 1) +
    facet_wrap(~ FDR_class, nrow = 4) +
    labs(
        x = "abs(Correlation Score)", y = "Density", fill = "Cell Type",
        color = "Cell Type"
    ) +
    theme_bw(base_size = 20)
pdf(file.path(plot_dir, "off_target_cor_stratified.pdf"), width = 10, height = 10)
print(p)
dev.off()

session_info()
