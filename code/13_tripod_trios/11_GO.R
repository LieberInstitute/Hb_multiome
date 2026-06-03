library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
library(Signac)
library(qs2)

cell_type_res = c('broad', 'fine')[
    as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
]

if (cell_type_res == 'broad') {
    cell_types = c(
        "Astrocyte", "Ependymal", "Excit.Thal", "Inhib.Thal",
        "Microglia", "Oligo", "OPC", "MHb", "LHb", "Inhib_LHb"
    )
} else {
    cell_types = c(
        'MHb.1', 'MHb.1.2', 'MHb.2', 'LHb.1.3.4', 'LHb.2.7', 'LHb.4',
        "Inhib_LHb_4.1", "Inhib_LHb_4.2", 'Excit.Thal', 'Inhib.Thal',
        'Astrocyte', 'Microglia', 'Oligo', 'OPC', 'Ependymal'
    )
}

seur_path = here(
    "processed-data", "13_tripod_trios", "02_tripod_preprocess",
    "preprocessed_objects_%s.qs2"
)
trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    sprintf("filtered_trios_%s.parquet", cell_type_res)
)
plot_dir = here("plots", "13_tripod_trios", "11_GO")
go_num_terms = 2

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

trio_df = read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    filter(is_unique) |>
    dplyr::select(cell_type, gene, TF) |>
    collect()

stopifnot(setequal(trio_df$cell_type, cell_types))

#   Run enrichGO per cell type
ego_df_list = list()
size_df_list = list()
for (this_cell_type in unique(trio_df$cell_type)) {
    this_trio_df = trio_df |>
        filter(cell_type == this_cell_type)
    gene_set = union(this_trio_df$gene, this_trio_df$TF)

    size_df_list[[this_cell_type]] = tibble(
        cell_type = this_cell_type,
        num_genes = length(gene_set)
    )
  
    #   Universe: all genes in the experiment
    universe_genes = rownames(
        qs_read(sprintf(seur_path, this_cell_type))$seur[['RNA']]
    )
    stopifnot(all(this_trio_df$gene %in% universe_genes))
    stopifnot(all(this_trio_df$TF %in% universe_genes))

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

message("Distribution of gene-set sizes across cell types:")
bind_rows(size_df_list) |>
    pull(num_genes) |>
    summary() |>
    print()

#   Custom dot plot by cell type
plot_df = bind_rows(ego_df_list) |>
    group_by(cell_type) |>
    slice_min(p.adjust, n = go_num_terms, with_ties = FALSE) |>
    ungroup() |>
    mutate(
        gene_ratio = Count / as.integer(str_extract(GeneRatio, "(?<=/)[0-9]+")),
        log_fdr    = -log10(p.adjust)
    )

#   Order GO terms by the first (highest-ranked) cell type they appear in
term_order = plot_df |>
    mutate(cell_type = factor(cell_type, levels = cell_types)) |>
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

pdf(file.path(plot_dir, "GO_trio_genes.pdf"), height = 5)
print(p)
dev.off()

session_info()
