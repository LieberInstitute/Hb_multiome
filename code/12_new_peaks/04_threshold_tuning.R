#   Try to tune FDR and correlation thresholds for linked peaks in a data-driven
#   way. The assumption here is that the best FDR and correlation thresholds will
#   yield the largest number of "real" results, while maintaining an acceptable
#   number of false positives. We assume that linked genes that are unique to a
#   given cell type are more likely than random to be markers for that cell type.
#   We'll tune thresholds to maximize the odds ratio for uniquely linked genes
#   being enriched for markers for the target cell type, provided a minimum
#   count of linked genes is met (to avoid noise in the OR metric). This metric
#   will be averaged across cell types

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
marker_path = here(
    'processed-data', '10_MAGMA', 'RNA', 'registration_banksy',
    'modeling_results', 'fine.rds'
)
plot_dir = here("plots", "12_new_peaks")
min_num_links = 10
FDR_thresholds = c(0.01, 0.05, 0.1, 0.15, 0.2, 1)
cor_thresholds = 0.05 * seq(0, 10)
marker_FDR = 0.1
max_markers = 100

################################################################################
#   Functions
################################################################################

is_unique_enough = function(target_cell_type, other_cell_type) {
    unique_enough = identical(target_cell_type, unique(other_cell_type)) ||
        (
            grepl('^[ML]Hb', target_cell_type) &&
            all(grepl('^[ML]Hb', other_cell_type))
        )
    return(unique_enough)
}

################################################################################
#   Import markers
################################################################################

marker_df = readRDS(marker_path)$enrichment |>
    as_tibble() |>
    pivot_longer(
        cols = matches('^(t_stat|p_value|fdr|logFC)_'),
        names_to = c('.value', 'cell_type'),
        names_pattern = '^(t_stat|p_value|fdr|logFC)_(.+)$'
    ) |>
    filter(fdr < marker_FDR, logFC > 0) |>
    group_by(cell_type) |>
    arrange(fdr) |>
    slice_head(n = max_markers) |>
    ungroup() |>
    select(cell_type, gene)

################################################################################
#   Collect odds ratio for enrichment of markers for each cell type
################################################################################

metric_df_list = list()
for (target_cell_type in cell_types) {
    result_df_list = list()
    for (other_cell_type in cell_types) {
        result_df_list[[other_cell_type]] = read_csv(
            sprintf(result_paths, target_cell_type, other_cell_type),
            show_col_types = FALSE
        ) |>
        select(peak, gene, score, FDR, other_cell_type)
    }
    result_df = bind_rows(result_df_list)

    for (FDR_threshold in FDR_thresholds) {
        for (cor_threshold in cor_thresholds) {
            linked_genes = result_df |>
                #   Filter all linked peaks by thresholds
                filter(FDR <= FDR_threshold, abs(score) >= cor_threshold) |>
                #   Take linked peaks now present only in the target cell type,
                #   or only among Hb cell types if the target is a Hb cell type
                group_by(peak, gene) |>
                filter(is_unique_enough(target_cell_type, other_cell_type)) |>
                ungroup() |>
                #   Then grab their unique genes
                pull(gene) |>
                unique()

            these_markers = marker_df |>
                filter(cell_type == target_cell_type) |>
                pull(gene)
            
            # Get all genes that could potentially be linked
            all_possible_genes = result_df |>
                pull(gene) |>
                unique()
            
            # Build 2x2 contingency table
            # Rows: is_linked (yes/no), Cols: is_marker (yes/no)
            in_linked_and_marker = sum(linked_genes %in% these_markers)
            in_linked_not_marker = length(linked_genes) - in_linked_and_marker
            not_linked_but_marker = sum(all_possible_genes %in% these_markers) - in_linked_and_marker
            not_linked_not_marker = length(all_possible_genes) - length(linked_genes) - not_linked_but_marker
            
            contingency_table = matrix(
                c(
                    in_linked_and_marker, in_linked_not_marker,
                    not_linked_but_marker, not_linked_not_marker
                ),
                dimnames = list(
                    c("Linked", "Not Linked"),
                    c("Marker", "Not Marker")
                ),
                nrow = 2, byrow = TRUE
            )
            
            fisher_result = fisher.test(contingency_table, alternative = "greater")
            
            metric_df_list[[length(metric_df_list) + 1]] = tibble(
                target_cell_type = target_cell_type,
                FDR_threshold = !!FDR_threshold,
                cor_threshold = !!cor_threshold,
                num_linked_genes = length(linked_genes),
                fisher_p = fisher_result$p.value,
                fisher_OR = fisher_result$estimate
            )
        }
    }
}
metric_df = bind_rows(metric_df_list)

################################################################################
#   Plot heatmaps
################################################################################

metric_df = metric_df |>
    mutate(
        FDR_threshold = factor(
            FDR_threshold, levels = sort(unique(FDR_threshold))
        ),
        fisher_OR = ifelse(
            num_linked_genes >= min_num_links, fisher_OR, NaN
        )
    )

p = ggplot(
        metric_df,
        aes(x = FDR_threshold, y = cor_threshold, fill = fisher_OR)
    ) +
    geom_tile() +
    scale_fill_viridis_c() +
    facet_wrap(~target_cell_type) +
    labs(
        x = "FDR Threshold",
        y = "Correlation Threshold",
        fill = "Enrichment OR",
        title = "Marker Enrichment Across Thresholds"
    ) +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
pdf(file.path(plot_dir, "threshold_heatmap_faceted.pdf"))
print(p)
dev.off()

p = metric_df |>
    group_by(FDR_threshold, cor_threshold) |>
    summarize(
        mean_prop_markers = ifelse(
            mean(is.na(prop_markers)) > 0.5,
            NaN,
            mean(prop_markers, na.rm = TRUE)
        )
    ) |>
    ggplot(aes(x = FDR_threshold, y = cor_threshold, fill = mean_prop_markers)) +
        geom_tile() +
        scale_fill_viridis_c() +
        labs(
            x = "FDR Threshold",
            y = "Correlation Threshold",
            fill = "Mean Proportion\nMarkers",
            title = "Mean Marker Enrichment Across Thresholds"
        ) +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
pdf(file.path(plot_dir, "threshold_heatmap_mean.pdf"))
print(p)
dev.off()

session_info()
