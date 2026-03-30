library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(duckplyr)
library(Polychrome)
library(cowplot)

seur_path = here(
    "processed-data", "06_peak_calling", "12_pseudobulk_MACS2",
    "Mid_pseudobulk.spearman.5e5.rds"
)
link_path = here(
    'processed-data', '12_new_peaks', '01_link_peaks', 'all_data.parquet'
)
plot_path = here(
    "plots", "12_new_peaks", "15_donor_specificity", "donor_specificity.pdf"
)
atac_assay = "ATAC_macs2_pseudo"
rna_assay = "RNA"
cor_thres = 0.3
FDR_thres = 0.1
num_expected_donors = 10
num_links = 20

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(plot_path), showWarnings = FALSE)

seur = readRDS(seur_path)
seur@meta.data$donor = str_extract(colnames(seur), 'S[0-9]{2}.*')
stopifnot(length(unique(seur@meta.data$donor)) == num_expected_donors)

donor_colors = palette36.colors(num_expected_donors)
names(donor_colors) = unique(seur@meta.data$donor)

link_df = read_parquet_duckdb(link_path) |>
    #   Note the lack of abs() is intentional here; we only want to visualize
    #   positive correlations in the heatmap
    filter(
        score > cor_thres, FDR < FDR_thres, target_cell_type == other_cell_type
    ) |>
    select(peak, gene, target_cell_type, FDR) |>
    #   Here I use a more relaxed definition of cell-type specificity. This is
    #   really one of two ways of doing this with the results that we have
    group_by(peak, gene) |>
    filter(n() == 1) |>
    ungroup() |>
    collect()

stopifnot(all(link_df$peak %in% rownames(seur[[atac_assay]])))
stopifnot(all(link_df$gene %in% rownames(seur[[rna_assay]])))

#   Sample linked peaks from top, bottom, and middle ranges of FDR by cell type
link_df_list = list()
for (this_percentile in c(0, 0.5)) {
    link_df_list[[as.character(this_percentile)]] = link_df |>
        group_by(target_cell_type) |>
        filter(
            FDR > quantile(FDR, this_percentile),
            n() > num_links * 4
        ) |>
        arrange(FDR) |>
        slice_head(n = num_links) |>
        ungroup() |>
        mutate(percentile = this_percentile)
}
link_df_list[['1']] = link_df |>
    group_by(target_cell_type) |>
    filter(n() > num_links * 4) |>
    arrange(desc(FDR)) |>
    slice_head(n = num_links) |>
    ungroup() |>
    mutate(percentile = 1)
link_sub_df = bind_rows(link_df_list)

count_df_list = list()
for (this_donor in unique(seur@meta.data$donor)) {
    subset_vec = (seur@meta.data$donor == this_donor)
    
    count_df_list[[length(count_df_list) + 1]] = tibble(
        cell_type = link_sub_df$target_cell_type,
        donor = this_donor,
        percentile = link_sub_df$percentile,
        peak = link_sub_df$peak,
        gene = link_sub_df$gene,
        peak_value = rowMeans(
            GetAssayData(
                seur, assay = atac_assay, layer = "data"
            )[link_sub_df$peak, subset_vec, drop = FALSE]
        ),
        gene_value = rowMeans(
            GetAssayData(
                seur, assay = rna_assay, layer = "data"
            )[link_sub_df$gene, subset_vec, drop = FALSE]
        )
    )
}
count_df = bind_rows(count_df_list) |>
    group_by(peak, gene) |>
    #   Z-scoring across donors helps show the influence of donor independent of
    #   the magnitiudes of peak or gene expression
    mutate(
        peak_value = (peak_value - mean(peak_value)) / sd(peak_value),
        gene_value = (gene_value - mean(gene_value)) / sd(gene_value)
    ) |>
    ungroup()

cell_types <- unique(count_df$cell_type)
percentiles <- unique(count_df$percentile)

# Generate one plot per cell_type x percentile
plots <- lapply(cell_types, function(ct) {
    lapply(percentiles, function(p) {
        df_sub <- count_df |>
            filter(cell_type == ct, percentile == p)
        
        ggplot(df_sub, aes(x = peak_value, y = gene_value, color = donor)) +
            geom_point(size = 0.5) +
            scale_color_manual(values = donor_colors) +
            labs(title = sprintf("%s | p=%.1f", ct, p), x = "Peak", y = "Gene") +
            theme_bw(base_size = 15) +
            theme(legend.position = "none")
    })
})

# Flatten to a single list, row-major (cell_type changes slowly)
plot_list <- unlist(plots, recursive = FALSE)

# Extract shared legend
legend <- get_legend(
    ggplot(count_df, aes(x = peak_value, y = gene_value, color = donor)) +
        geom_point() +
        scale_color_manual(values = donor_colors) +
        theme_bw(base_size = 10) +
        guides(color = guide_legend(override.aes = list(size = 2)))
)

grid <- plot_grid(plotlist = plot_list, ncol = length(percentiles))
final <- plot_grid(grid, legend, rel_widths = c(1, 0.15))

pdf(plot_path, width = 10, height = 3 * length(cell_types))
print(final)
dev.off()
