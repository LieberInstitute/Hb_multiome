#   Try to tune FDR and correlation thresholds for linked peaks in a data-driven
#   way. The assumption here is that the best FDR and correlation thresholds will
#   yield the largest number of "real" results, while maintaining an acceptable
#   number of false positives. We assume that linked genes that are unique to a
#   given cell type are more likely than random to be markers for that cell type.
#   We'll tune thresholds to maximize significance for uniquely linked genes
#   being enriched for markers for the target cell type

library(here)
library(tidyverse)
library(sessioninfo)
library(duckplyr)

cell_types = c(
    'Astrocyte', 'Endo', 'Excit.Thal', 'Inhib.Thal', 'LHb.1', 'LHb.1.3',
    'LHb.1.3.4', 'LHb.2.7', 'LHb.4', 'LHb.7', 'MHb.1', 'MHb.1.2', 'MHb.2',
    'MHb.3', 'Microglia', 'Oligo', 'OPC', 'Thal'
)
link_path = here(
    'processed-data', '12_new_peaks', '01_link_peaks', 'all_data.parquet'
)
marker_path = here(
    'processed-data', '10_MAGMA', 'RNA', 'registration_banksy',
    'modeling_results', 'fine.rds'
)
out_path = here(
    'processed-data', '12_new_peaks', '04_threshold_tuning',
    'marker_enrichment_metrics.csv'
)
plot_dir = here("plots", "12_new_peaks", "04_threshold_tuning")
FDR_thresholds = c(0.01, 0.05, 0.1, 0.15, 0.2, 1)
cor_thresholds = 0.05 * seq(0, 10)
marker_FDR = 0.1
max_markers = 100

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(out_path), showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

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

all_possible_genes = read_parquet_duckdb(link_path, prudence = 'stingy') |>
    distinct(gene) |>
    collect() |>
    pull(gene)

metric_df_list = list()
for (FDR_threshold in FDR_thresholds) {
    for (cor_threshold in cor_thresholds) {
        link_unique_df = read_parquet_duckdb(link_path, prudence = 'lavish') |>
            filter(score > abs(cor_threshold), FDR < FDR_threshold) |>
            distinct(peak, gene, other_cell_type) |>
            group_by(peak, gene) |>
            filter(
                (n() == 1) | all(grepl('^[ML]Hb', unique(other_cell_type)))
            ) |>
            ungroup() |>
            dplyr::rename(cell_type = other_cell_type) |>
            collect()

        for (this_cell_type in cell_types) {
            linked_genes = link_unique_df |>
                filter(cell_type == this_cell_type) |>
                pull(gene) |>
                unique()
            
            these_markers = marker_df |>
                filter(cell_type == this_cell_type) |>
                pull(gene)
          
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
                cell_type = this_cell_type,
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

#   Save in case we want to plot interactively (the above takes a long time to
#   compute)
write_csv(metric_df, out_path)

################################################################################
#   Plot heatmaps
################################################################################

metric_df = metric_df |>
    mutate(
        FDR_threshold = factor(
            FDR_threshold, levels = sort(unique(FDR_threshold))
        ),
        fisher_log10p = -log10(fisher_p),
        fisher_logOR = replace_na(log(fisher_OR + 1), 0)
    )

for (this_metric in c("fisher_log10p", "fisher_logOR")) {
    if (this_metric == "fisher_log10p") {
        fill_title = "Enrichment\n-log10(p-value)"
    } else if (this_metric == "fisher_logOR") {
        fill_title = "Enrichment\nlog(OR + 1)"
    }
  
    p = ggplot(
            metric_df,
            aes(x = FDR_threshold, y = cor_threshold, fill = !!sym(this_metric))
        ) +
        geom_tile() +
        scale_fill_viridis_c() +
        facet_wrap(~cell_type) +
        labs(
            x = "FDR Threshold",
            y = "Correlation Threshold",
            fill = fill_title,
            title = "Marker Enrichment Across Thresholds"
        ) +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    pdf(file.path(plot_dir, sprintf("faceted_%s.pdf", this_metric)))
    print(p)
    dev.off()

    p = metric_df |>
        group_by(FDR_threshold, cor_threshold) |>
        summarize(mean_metric = mean(!!sym(this_metric))) |>
        ggplot(aes(x = FDR_threshold, y = cor_threshold, fill = mean_metric)) +
            geom_tile() +
            scale_fill_viridis_c() +
            labs(
                x = "FDR Threshold",
                y = "Correlation Threshold",
                fill = paste("Mean", fill_title)
            ) +
            theme_bw(base_size = 15) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    pdf(file.path(plot_dir, sprintf("mean_%s.pdf", this_metric)))
    print(p)
    dev.off()
}

session_info()
  