library(here)
library(tidyverse)
library(sessioninfo)
library(clusterProfiler)
library(duckplyr)
library(org.Hs.eg.db)

peak_path = here(
    'processed-data', '12_new_peaks', '01_link_peaks', 'all_data.csv.gz'
)
plot_dir = here("plots", "12_new_peaks", "11_peak_go")
fdr_cutoff_peak = 0.1
fdr_cutoff_go = 0.05
cor_cutoff = 0.3

dir.create(plot_dir, showWarnings = FALSE)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))

################################################################################
#   Functions
################################################################################

do_go = function(gene_list, universe, plot_path, num_terms = 2) {
    go_obj = compareCluster(
        gene_list, fun = "enrichGO", universe = universe,
        OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE,
        keyType = "SYMBOL"
    )

    if(!is.null(go_obj)) {
        go_obj@compareClusterResult = go_obj@compareClusterResult |>
            filter(p.adjust < fdr_cutoff_go)
        
        if(nrow(go_obj@compareClusterResult) > 0) {
            pdf(plot_path)
            print(dotplot(go_obj, showCategory = num_terms))
            dev.off()
        }
    }
}

################################################################################
#   Main
################################################################################

peak_df = read_csv_duckdb(peak_path, prudence = 'lavish')
background_universe = unique(peak_df$gene)

peak_df = peak_df |>
    filter(FDR < fdr_cutoff_peak, abs(score) >= cor_cutoff) |>
    group_by(peak, gene, target_cell_type) |>
    mutate(is_unique = unique(other_cell_type) == target_cell_type) |>
    ungroup() |>
    filter(target_cell_type == other_cell_type)

for (cell_type in unique(peak_df$target_cell_type)) {
    this_peak_df = peak_df |>
        filter(target_cell_type == cell_type)

    gene_list = list()
    gene_list[['shared_positive']] = this_peak_df |>
        filter(score > 0, !is_unique) |>
        pull(gene) |>
        unique()
    gene_list[['unique_positive']] = this_peak_df |>
        filter(score > 0, is_unique) |>
        pull(gene) |>
        unique()
    gene_list[['shared_negative']] = this_peak_df |>
        filter(score < 0, !is_unique) |>
        pull(gene) |>
        unique()
    gene_list[['unique_negative']] = this_peak_df |>
        filter(score < 0, is_unique) |>
        pull(gene) |>
        unique()

    do_go(
        gene_list, background_universe,
        file.path(plot_dir, sprintf("%s_GO.pdf", cell_type))
    )
}

session_info()
