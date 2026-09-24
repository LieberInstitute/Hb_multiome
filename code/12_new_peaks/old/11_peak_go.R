library(here)
library(tidyverse)
library(sessioninfo)
library(clusterProfiler)
library(duckplyr)
library(org.Hs.eg.db)

link_all_path = here(
    'processed-data', '12_new_peaks', '01_link_peaks', 'all_data.parquet'
)
link_filtered_path = here(
    'processed-data', '12_new_peaks', '01_link_peaks', 'filtered_data.parquet'
)
plot_dir = here("plots", "12_new_peaks", "11_peak_go")
fdr_cutoff_peak = 0.1
fdr_cutoff_go = 0.1
cor_cutoff = 0.3

dir.create(plot_dir, showWarnings = FALSE)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

################################################################################
#   Functions
################################################################################

do_go = function(gene_list, universe, plot_path, num_terms = 3) {
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
            pdf(plot_path, width = 5)
            print(dotplot(go_obj, showCategory = num_terms))
            dev.off()
        }
    }
}

################################################################################
#   Main
################################################################################

link_df = read_parquet_duckdb(link_all_path, prudence = 'lavish')
background_universe = unique(link_df$gene)

#   I tried splitting by unique vs. shared links, but almost none were shared.
#   Here just use all pairs together
link_df = read_parquet_duckdb(link_filtered_path, prudence = 'stingy') |>
    collect()

for (this_cell_type in unique(link_df$cell_type)) {
    this_link_df = link_df |>
        filter(cell_type == this_cell_type)

    gene_list = list()
    gene_list[['positive']] = this_link_df |>
        filter(score > 0) |>
        pull(gene) |>
        unique()
    gene_list[['negative']] = this_link_df |>
        filter(score < 0) |>
        pull(gene) |>
        unique()

    do_go(
        gene_list, background_universe,
        file.path(plot_dir, sprintf("%s_GO.pdf", this_cell_type))
    )
}

session_info()
