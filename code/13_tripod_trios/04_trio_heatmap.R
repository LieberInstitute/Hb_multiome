library(sessioninfo)
library(tidyverse)
library(here)
library(ComplexHeatmap)
library(circlize)
library(duckplyr)
library(Seurat)
library(Signac)
library(qs2)
library(GenomicRanges)

trio_paths = here(
    "processed-data", "12_new_peaks", "10_tripod_trios", "trios_%s.csv.gz"
)
# trio_paths = here(
#     "processed-data", "13_tripod_trios", "03_tripod_trios", "trios_%s.parquet"
# )
link_filtered_path = here(
    "processed-data", "12_new_peaks", "01_link_peaks", "old",
    "filtered_data.parquet"
)
# link_filtered_path = here(
#     "processed-data", "12_new_peaks", "01_link_peaks", "filtered_data.parquet"
# )
link_all_path = here(
    "processed-data", "12_new_peaks", "01_link_peaks", "old",
    "all_data.parquet"
)
# link_all_path = here(
#     "processed-data", "12_new_peaks", "01_link_peaks",
#     "all_data.parquet"
# )
seur_path = here(
    "processed-data", "13_tripod_trios", "01_chromVAR", "seur.qs2"
)
# seur_path = here(
#     'processed-data', '11_link_prep', '02_rebuild_atac_assay',
#     'cell_level_seur.qs2'
# )
plot_path = here(
    "plots", "13_tripod_trios", "04_trio_heatmap", "heatmap.pdf"
)
cell_type1 = "MHb.2"
cell_type2 = "LHb.2.7"
cor_thres = 0.3
FDR_thres_trio = 0.1
FDR_thres_link = 0.1

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(plot_path), showWarnings = FALSE)

link_filtered_df = read_parquet_duckdb(
        link_filtered_path, prudence = 'stingy'
    ) |>
    collect()

trio_df_list = list()
for (cell_type in c(cell_type1, cell_type2)) {
    trio_df_list[[cell_type]] = read_csv_duckdb(
            sprintf(trio_paths, cell_type), prudence = 'lavish'
        ) |>
        filter(
            condition_on == 'Yj', stringency_level == 2, adj < FDR_thres_trio
        ) |>
        select(peak, gene, TF, adj) |>
        mutate(cell_type = cell_type)
}
trio_df = bind_rows(trio_df_list) |>
    collect()

final_df = trio_df |>
    inner_join(
        link_filtered_df |>
            select(peak, gene, cell_type),
        by = c("peak", "gene", "cell_type")
    )

seur = qs_read(seur_path)
Idents(seur) = seur@meta.data$refined_mid_cluster

#   Really just checks we're importing the right data, as this certainly
#   should be true
stopifnot(all(final_df$peak %in% rownames(seur[['ATAC']])))
stopifnot(all(final_df$gene %in% rownames(seur[['RNA']])))

#   Import links for the cell type of interest
a = here(
    "processed-data", "12_new_peaks", "01_link_peaks", "old",
    "LHb.2.7_LHb.2.7.csv.gz"
)
Links(seur[['ATAC']]) = read_csv_duckdb(a, prudence = "stingy") |>
    filter(FDR < 0.1, abs(score) > 0.3) |>
    select(
        seqnames, start, end, width, strand, score, gene, peak, zscore, pvalue
    ) |>
    as_tibble() |>
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

this_gene = final_df$gene[1]
this_peak = final_df$peak[1]
this_TF = final_df$TF[1]

#   Determine an appropriate region to plot by forming a small window
#   around both the gene and link
gene_coords = LookupGeneCoords(seur, gene = this_gene)
all_positions = c(
    start(gene_coords), end(gene_coords),
    str_extract(this_peak, '-([0-9]+)-', group = 1) |> as.numeric(),
    str_extract(this_peak, '([0-9]+)$', group = 1) |> as.numeric()
)
buffer = as.integer(0.1 * (max(all_positions) - min(all_positions)))
if (max(all_positions) == end(gene_coords)) {
    extend_downstream = buffer
    extend_upstream = start(gene_coords) - min(all_positions) + buffer
} else {
    extend_upstream = buffer
    extend_downstream = max(all_positions) - end(gene_coords) + buffer
}

p = CoveragePlot(
        seur, region = final_df$gene[4], extend.upstream = extend_upstream,
        extend.downstream = extend_downstream, links = final_df$gene[4],
        region.highlight = StringToGRanges(final_df$peak[4])
    ) & 
    theme(text = element_text(size = 14))
p = p[[1]] +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


