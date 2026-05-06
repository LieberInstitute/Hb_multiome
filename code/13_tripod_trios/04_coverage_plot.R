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
    "processed-data", "13_tripod_trios", "03_tripod_trios", "trios_%s.parquet"
)
link_filtered_path = here(
    'processed-data', '12_new_peaks', '15_donor_specificity',
    'metacell_filtered_unbiased.parquet'
)
link_cell_type_path = here(
    "processed-data", "12_new_peaks", "14_metacell_link_peaks", "%s.parquet"
)
preprocessed_path = here(
    "processed-data", "13_tripod_trios", "02_tripod_preprocess",
    "preprocessed_objects_%s.qs2"
)
seur_path = here(
    "processed-data", "13_tripod_trios", "01_chromVAR", "seur.qs2"
)
plot_path = here(
    "plots", "13_tripod_trios", "04_coverage_plot", "coverage_plot.pdf"
)
cell_type1 = "MHb.2"
cell_type2 = "LHb.2.7"
FDR_thres_trio = 0.05

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(plot_path), showWarnings = FALSE, recursive = TRUE)

################################################################################
#   Import and filter trios, overlapping with existing links
################################################################################

link_filtered_df = read_parquet_duckdb(
        link_filtered_path, prudence = 'lavish'
    ) |>
    filter(score > 0) |>
    collect()

trio_df_list = list()
for (cell_type in c(cell_type1, cell_type2)) {
    trio_df_list[[cell_type]] = read_parquet_duckdb(
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
    inner_join(link_filtered_df, by = c("peak", "gene", "cell_type"))

this_gene = final_df$gene[1]
this_peak = final_df$peak[1]
this_TF = final_df$TF[1]
this_cell_type = final_df$cell_type[1]

################################################################################
#   Import Seurat object and attach links for the relevant cell type
################################################################################

seur = qs_read(seur_path)
Idents(seur) = seur@meta.data$refined_mid_cluster

#   Really just checks we're importing the right data, as this certainly
#   should be true
stopifnot(all(final_df$peak %in% rownames(seur[['ATAC']])))
stopifnot(all(final_df$gene %in% rownames(seur[['RNA']])))

#   Import links for the cell type of interest
Links(seur[['ATAC']]) = sprintf(link_cell_type_path, this_cell_type) |>
    read_parquet_duckdb(prudence = "stingy") |>
    filter(FDR < FDR_thres_link, abs(score) > cor_thres) |>
    select(
        seqnames, start, end, width, strand, score, gene, peak, zscore, pvalue
    ) |>
    collect() |>
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

################################################################################
#   Determine plotting region and parameters for the CoveragePlot
################################################################################

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

#   The actual GRanges used in the coverage plot. For some reason providing
#   this directly causes issues, but we'll need it for the motif track
real_window = sprintf(
        '%s-%d-%d',
        seqnames(gene_coords),
        min(all_positions) - extend_upstream,
        max(all_positions) + extend_downstream
    ) |>
    StringToGRanges()

################################################################################
#   Add a motif track manually
################################################################################

#   In TRIPOD, a TF is associated with exactly one motif. Grab the input objects
#   from TRIPOD to find that motif, and which peaks contain it. The motif track
#   will just show binary presence of the motif in measured peaks within the
#   plotting region

input_objs = qs_read(sprintf(preprocessed_path, this_cell_type))

motif_tf_map = input_objs$tripod_seur$motifxTF
stopifnot(this_TF %in% motif_tf_map[, 'TF'])
this_motif = motif_tf_map[motif_tf_map[, 'TF'] == this_TF, 'motif']

motif_present = input_objs$tripod_seur$peakxmotif[, this_motif, drop = TRUE]
motif_peaks = names(motif_present)[motif_present]
stopifnot(all(motif_peaks %in% rownames(seur[['ATAC']])))

motif_gr = granges(seur[['ATAC']])[rownames(seur[['ATAC']]) %in% motif_peaks, ]
hits = findOverlaps(motif_gr, real_window)
motif_gr = motif_gr[queryHits(hits)]

motif_track = tibble(start = start(motif_gr), end = end(motif_gr)) |>
    ggplot() +
        geom_rect(
            aes(xmin = start, xmax = end, ymin = 0, ymax = 1),
            fill = "forestgreen"
        ) +
        scale_x_continuous(
            limits = c(start(real_window), end(real_window))#,
            # expand = c(0, 0)
        ) +
        theme_void()

################################################################################
#   Construct the final coverage plot
################################################################################

p = (
        CoveragePlot(
            seur, region = this_gene, extend.upstream = extend_upstream,
            extend.downstream = extend_downstream, links = this_gene,
            region.highlight = StringToGRanges(this_peak)
        ) & 
        theme(text = element_text(size = 14))
    ) /
    motif_track
p = p +
    patchwork::plot_layout(heights = c(10, 0.7))
# p = p[[1]] +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(plot_path)
print(p)
dev.off()

session_info()
