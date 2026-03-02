#   Try to tune FDR and correlation thresholds for linked peaks in a data-driven
#   way. The assumption here is that the best FDR and correlation thresholds will
#   yield the largest number of "real" results, while maintaining an acceptable
#   number of false positives. We assume that linked genes that are unique to a
#   given cell type are more likely than random to be markers for that cell type.
#   We'll tune thresholds to maximize the proportion of uniquely linked genes
#   that are markers for the target cell type, provided a minimum count of
#   linked genes is met (to avoid noise in the proportion metric). This metric
#   will be averaged across cell types

library(here)
library(tidyverse)
library(sessioninfo)
library(rtracklayer)

cell_types = c(
    'Astrocyte', 'Endo', 'Excit.Thal', 'Inhib.Thal', 'LHb.1', 'LHb.1.3',
    'LHb.1.3.4', 'LHb.2.7', 'LHb.4', 'LHb.7', 'MHb.1', 'MHb.1.2', 'MHb.2',
    'MHb.3', 'Microglia', 'Oligo', 'OPC', 'Thal'
)
result_paths = here(
    'processed-data', '12_new_peaks', '01_link_peaks', '%s_%s.csv.gz'
)
gtf_path = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'
marker_path = here("processed-data", "10_MAGMA", "RNA", "gene_sets", "mid.tsv")
plot_path = here("plots", "12_new_peaks", "threshold_heatmap.pdf")
min_num_links = 10
FDR_thresholds = c(0.01, 0.05, 0.1, 0.15, 0.2, 1)
cor_thresholds = 0.05 * seq(0, 10)

################################################################################
#   Import markers and convert to gene symbols
################################################################################

gtf = import(gtf_path) |>
    as.data.frame() |>
    as_tibble() |>
    filter(type == 'gene') |>
    select(gene_id, gene_name)

marker_df = read_tsv(marker_path, show_col_types = FALSE) |>
    left_join(gtf, by = 'gene_id')

message(
    sprintf(
        "Proportion of markers with NA gene_name: %.1f%%",
        100 * mean(is.na(marker_df$gene_name))
    )
)

marker_df = marker_df |>
    filter(!is.na(gene_name)) |>
    select(set_id, gene_name)

################################################################################
#   Collect proportion of markers for each cell type
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
                #   Take linked peaks now present only in the target cell type
                group_by(peak, gene) |>
                filter(identical(target_cell_type, unique(other_cell_type))) |>
                ungroup() |>
                #   Then grab their unique genes
                pull(gene) |>
                unique()

            these_markers = marker_df |>
                filter(set_id == target_cell_type) |>
                pull(gene_name)

            metric_df_list[[length(metric_df_list) + 1]] = tibble(
                target_cell_type = target_cell_type,
                FDR_threshold = !!FDR_threshold,
                cor_threshold = !!cor_threshold,
                num_linked_genes = length(linked_genes),
                prop_markers = mean(linked_genes %in% these_markers)
            )
        }
    }
}
metric_df = bind_rows(metric_df_list)

################################################################################
#   Collect proportion of markers for each cell type
################################################################################

p = metric_df |>
    mutate(
        FDR_threshold = factor(
            FDR_threshold, levels = sort(unique(FDR_threshold))
        )
    ) |>
    ggplot(aes(x = FDR_threshold, y = cor_threshold, fill = prop_markers)) +
        geom_tile() +
        scale_fill_viridis_c() +
        facet_wrap(~target_cell_type) +
        labs(
            x = "FDR Threshold",
            y = "Correlation Threshold",
            fill = "Proportion\nMarkers",
            title = "Marker Enrichment Across Thresholds"
        ) +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
pdf(plot_path)
print(p)
dev.off()

session_info()
