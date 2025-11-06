library(tidyverse)
library(here)
library(getopt)
library(org.Hs.eg.db)
library(clusterProfiler)

pair_path = here(
    'processed-data', '06_peak_calling', '21_overlaping_FDRscores_TopHeatmap',
    'overlaps_linkPeak_DARs_classified_thr_CC0.3_thr_DAR0.1.csv'
)
universe_dir = here(
    'processed-data', '06_peak_calling',
    '13_pseudobulk_LinkPeaks_MACS2_split_ct', 'links_ct_merged', 'seurats_rds'
)

pairs_df = read_csv(pair_path, show_col_types = FALSE) |>
    filter(type_classification == 'classification-3') |>
    select(cell_type, gene_id, gene_name, peak_id_links, category)

#   Actually, this approach doesn't work, as they appear to all be the same gene
#   sets
universe_df_list = list()
for (cell_type in unique(pairs_df$cell_type)) {
    rds_path = file.path(
        universe_dir,
        sprintf('Mid_%s_pseudobulk_seurat_subset.rds', cell_type)
    )
    universe_df_list[[cell_type]] = tibble(
        cell_type = cell_type,
        gene_name = rownames(readRDS(rds_path))
    )
}
universe_df = bind_rows(universe_df_list)
