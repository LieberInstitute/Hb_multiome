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
library(patchwork)

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
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

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(plot_path), showWarnings = FALSE, recursive = TRUE)

################################################################################
#   Import and filter trios, overlapping significant level 1 and 2 
################################################################################

trio_df = read_parquet_duckdb(trio_path, prudence = "stingy") |>
    filter(is_intersect, stringency_level == 1, cell_type == cell_type1) |>
    collect()

link_gr = read_parquet_duckdb(trio_path, prudence = "stingy") |>
    filter(stringency_level == 1, cell_type == cell_type1) |>
    collect() |>
    dplyr::rename(score = coef, pvalue = adj) |>
    separate(peak, into = c("seqnames", "start", "end"), sep = "-", convert = TRUE) |>
    mutate(width = end - start, strand = "*") |>
    select(seqnames, start, end, width, strand, score, gene, pvalue) |>
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

this_gene = trio_df$gene[1]
this_peak = trio_df$peak[1]
this_TF = trio_df$TF[1]
this_cell_type = trio_df$cell_type[1]

################################################################################
#   Import Seurat object and attach links for the relevant cell type
################################################################################

seur = qs_read(seur_path)
Idents(seur) = case_when(
    seur@meta.data$refined_mid_cluster == cell_type1 ~ cell_type1,
    grepl('Hb', seur@meta.data$refined_mid_cluster) ~ 'Other habenula',
    grepl('Thal', seur@meta.data$refined_mid_cluster) ~ 'Thalamus',
    TRUE ~ 'Glia'
)

#   Really just checks we're importing the right data, as this certainly
#   should be true
stopifnot(all(trio_df$peak %in% rownames(seur[['ATAC']])))
stopifnot(all(trio_df$gene %in% rownames(seur[['RNA']])))

#   Import links for the cell type of interest
Links(seur[['ATAC']]) = link_gr

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
#   Construct the base CoveragePlot and extract its x-axis limits
################################################################################

# Build the base CoveragePlot without built-in links first, so we can extract
# the actual rendered x-axis range and use it for all custom tracks below.
cp_no_links = CoveragePlot(
    seur, region = this_gene, extend.upstream = extend_upstream,
    extend.downstream = extend_downstream,
    region.highlight = StringToGRanges(this_peak),
    links = FALSE, features = this_gene
)

# Extract the actual x-axis range from the first panel of the CoveragePlot.
# Using this for all custom tracks ensures they share the same coordinate space.
cp_built = ggplot_build(cp_no_links$patches$plots[[1]])
xlim = cp_built$layout$get_scales(1)$x$range$range

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

# Use xlim (from the CoveragePlot) rather than real_window so this track
# aligns with the other panels
motif_track = tibble(start = start(motif_gr), end = end(motif_gr)) |>
    ggplot() +
        geom_rect(
            aes(xmin = start, xmax = end, ymin = 0, ymax = 1),
            fill = "forestgreen"
        ) +
        coord_cartesian(xlim = xlim) +
        theme_void()

################################################################################
#   Build the peak -> gene link track
################################################################################

# Filter to only the specific peak-gene pair of interest (this_peak -> this_gene).
gene_tss = start(gene_coords)
links_to_gene = as.data.frame(link_gr) |>
    mutate(
        peak_name = paste(seqnames, start, end, sep = "-"),
        peak_mid = (start + end) / 2,
        gene_tss = gene_tss,
        y = 0, yend = 0
    ) |>
    filter(gene == this_gene, peak_name == this_peak)

# Use coord_cartesian instead of scale_x_continuous(limits = ...) so the arc
# of geom_curve is not clipped at the data level.
link_track_aligned = ggplot(links_to_gene, aes(x = peak_mid, xend = gene_tss, y = y, yend = yend)) +
  geom_curve(
    curvature = 0.2, ncp = 100,
    arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
    linewidth = 1, color = "steelblue"
  ) +
  coord_cartesian(xlim = xlim, clip = "off") +
  theme_void()

# Signac's CombineTracks internally places the expression violin to the right
# of the coverage tracks with widths = c(10, length(features)). Since our
# custom tracks are single panels, patchwork would stretch them to the full
# combined width. Adding a plot_spacer() with the same 10:1 ratio ensures
# the coverage portion aligns across all rows.
# this_gene is the features argument passed to CoveragePlot
cp_col_widths = c(10, length(this_gene))

link_track_row = (link_track_aligned | plot_spacer()) +
    plot_layout(widths = cp_col_widths)
motif_track_row = (motif_track | plot_spacer()) +
    plot_layout(widths = cp_col_widths)

# Combine with CoveragePlot and motif track, aligning axes
p = (cp_no_links & theme(text = element_text(size = 14))) /
    link_track_row /
    motif_track_row +
    plot_layout(heights = c(10, 1.5, 0.7), guides = "collect")

pdf(plot_path)
print(p)
dev.off()

session_info()
