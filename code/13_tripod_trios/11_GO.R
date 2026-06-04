library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)
library(Signac)
library(qs2)
library(rtracklayer)

cell_type_res = c('broad', 'fine')[
    as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
]

if (cell_type_res == 'broad') {
    cell_types = c("MHb", "LHb", "Inhib_LHb")
    go_num_terms = 4
} else {
    cell_types = c(
        'MHb.1', 'MHb.1.2', 'MHb.2', 'LHb.1.3.4', 'LHb.2.7', 'LHb.4',
        "Inhib_LHb_4.1", "Inhib_LHb_4.2", 'Excit.Thal', 'Inhib.Thal',
        'Astrocyte', 'Microglia', 'Oligo', 'OPC', 'Ependymal'
    )
    go_num_terms = 2
}

seur_path = here(
    "processed-data", "13_tripod_trios", "02_tripod_preprocess",
    "preprocessed_objects_%s.qs2"
)
trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    sprintf("filtered_trios_%s.parquet", cell_type_res)
)
out_path = here(
    "processed-data", "13_tripod_trios", "11_GO",
    sprintf("gene_sets_%s.tsv", cell_type_res)
)
gtf_path = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'
plot_dir = here("plots", "13_tripod_trios", "11_GO")
max_genes = 100

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

trio_df = read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    filter(is_unique) |>
    dplyr::select(cell_type, gene, TF, adj, coef, stringency_level) |>
    collect()

stopifnot(setequal(trio_df$cell_type, cell_types))

#   Run enrichGO per cell type
ego_df_list = list()
size_df_list = list()
gene_set_df_list = list()
for (this_cell_type in unique(trio_df$cell_type)) {
    gene_set = trio_df |>
        filter(cell_type == this_cell_type) |>
        pivot_longer(
            cols = c(gene, TF), names_to = "gene_type", values_to = "gene"
        ) |>
        #   If genes are duplicated, prioritize significance first then
        #   stringency level
        group_by(gene) |>
        arrange(adj, stringency_level) |>
        slice_head(n = 1) |>
        #   Then evenly sample genes and TFs among both stringency levels,
        #   prioritizing significance then effect size
        group_by(stringency_level, gene_type) |>
        arrange(adj, desc(coef)) |>
        slice_head(n = as.integer(max_genes / 4)) |>
        pull(gene)

    size_df_list[[this_cell_type]] = tibble(
        cell_type = this_cell_type,
        num_genes = length(gene_set)
    )
  
    gene_set_df_list[[this_cell_type]] = tibble(
        cell_type = this_cell_type,
        gene = gene_set
    )
  
    #   Universe: all genes expressed in the cell type
    universe_genes = rownames(
        qs_read(sprintf(seur_path, this_cell_type))$seur[['RNA']]
    )
    stopifnot(all(gene_set %in% universe_genes))

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

gtf = import(gtf_path) |>
    as.data.frame() |>
    as_tibble() |>
    filter(type == "gene") |>
    dplyr::rename(gene = gene_name) |>
    dplyr::select(gene_id, gene)
    
#   Will need gene sets (as ENSEMBL ID) for MAGMA
gene_set_df = bind_rows(gene_set_df_list) |>
    left_join(gtf, by = "gene") |>
    dplyr::rename(set_id = cell_type) |>
    dplyr::select(set_id, gene_id)

message(
    sprintf("Dropping %d genes not in the GTF", sum(is.na(gene_set_df$gene_id)))
)
gene_set_df |>
    filter(!is.na(gene_id)) |>
    write_tsv(out_path)

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
        cell_type   = factor(cell_type, levels = cell_types),
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

pdf(
    file.path(plot_dir, sprintf("GO_trio_genes_%s.pdf", cell_type_res)),
    width = 3 + length(cell_types) / 3, height = 5
)
print(p)
dev.off()

session_info()
