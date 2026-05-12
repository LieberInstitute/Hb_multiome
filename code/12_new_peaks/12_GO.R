library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

link_path = here(
    'processed-data', '12_new_peaks', '05_donor_specificity',
    'metacell_filtered_unbiased.parquet'
)
plot_dir = here("plots", "12_new_peaks", "12_GO")
cell_type_levels = c(
    'MHb.1', 'MHb.1.2', 'MHb.2', 'MHb.3', 'LHb.1.3.4', 'LHb.2.7', 'LHb.4',
    'Excit.Thal', 'Inhib.Thal', 'Astrocyte', 'Endo', 'Microglia', 'Oligo',
    'OPC', 'Ependymal'
)
go_num_terms = 2

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

link_df = read_parquet_duckdb(link_path, prudence = 'stingy') |>
    dplyr::select(gene, cell_type) |>
    collect()

#   Universe: all unique genes across all cell types
universe_genes = unique(link_df$gene)

#   Run enrichGO per cell type
ego_df_list = list()
for (this_cell_type in unique(link_df$cell_type)) {
    gene_set = link_df |>
        filter(cell_type == this_cell_type) |>
        pull(gene) |>
        unique()

    ego = enrichGO(
        gene          = gene_set,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        universe      = universe_genes,
        pAdjustMethod = "BH",
        pvalueCutoff  = 1,
        qvalueCutoff  = 0.05
    )
    ego_df_list[[this_cell_type]] = ego@result |>
        as_tibble() |>
        mutate(cell_type = this_cell_type)
}

#   Custom dot plot by cell type
plot_df = bind_rows(ego_df_list) |>
    group_by(cell_type) |>
    slice_min(p.adjust, n = go_num_terms, with_ties = FALSE) |>
    ungroup() |>
    mutate(
        gene_ratio = Count / as.integer(str_extract(GeneRatio, "(?<=/)[0-9]+")),
        log_fdr    = -log10(p.adjust)
    )

#   Validate all cell types in results are covered by cell_type_levels
stopifnot(all(unique(plot_df$cell_type) %in% cell_type_levels))

#   Order x-axis by cell_type_levels (keeping only those present)
ct_order = cell_type_levels[cell_type_levels %in% unique(plot_df$cell_type)]

#   Order GO terms by the first (highest-ranked) cell type they appear in
term_order = plot_df |>
    mutate(cell_type = factor(cell_type, levels = ct_order)) |>
    group_by(Description) |>
    slice_min(cell_type, n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(cell_type) |>
    pull(Description)

p = plot_df |>
    mutate(
        cell_type   = factor(cell_type, levels = ct_order),
        Description = factor(Description, levels = term_order)
    ) |>
    ggplot(
        aes(x = cell_type, y = Description, color = log_fdr, size = gene_ratio)
    ) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(
        x = "Cell Type", y = "GO Term", color = "-log10(FDR)",
        size = "Gene Ratio"
    )

pdf(file.path(plot_dir, "GO_linked_genes.pdf"), height = 5)
print(p)
dev.off()

session_info()
