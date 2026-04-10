#   How much are overall peak-gene correlations measured by linked peaks
#   influenced by a handful of donors? Produce various visualizations and
#   metrics to quantify this bias

library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(duckplyr)
library(Polychrome)
library(cowplot)
library(qs2)

dataset = c("pb", "metacell")[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

if (dataset == "pb") {
    result_paths = here(
        'processed-data', '12_new_peaks', '01_link_peaks', '%s.parquet'
    )
} else {
    result_paths = here(
        "processed-data", "12_new_peaks", "14_metacell_link_peaks", "%s.parquet"
    )
}

seur_path = here(
    'processed-data', '11_link_prep', '03_pseudobulk', 'pb_seur.qs2'
)
link_path = here(
    'processed-data', '12_new_peaks', '14_metacell_link_peaks',
    sprintf('%s_filtered_data.parquet', dataset)
)
plot_dir = here("plots", "12_new_peaks", "15_donor_specificity")
atac_assay = "ATAC"
rna_assay = "RNA"
num_expected_donors = 10
num_links = 20

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

seur = qs_read(seur_path)
seur@meta.data$donor = str_extract(colnames(seur), 'S[0-9]{2}.*')
stopifnot(length(unique(seur@meta.data$donor)) == num_expected_donors)

donor_colors = palette36.colors(num_expected_donors)
names(donor_colors) = unique(seur@meta.data$donor)
link_colors = palette36.colors(num_links)

link_df = read_parquet_duckdb(link_path) |>
    #   Cell-type specific, positively correlated links only
    filter(!is_shared, score > 0) |>
    collect()

stopifnot(all(link_df$peak %in% rownames(seur[[atac_assay]])))
stopifnot(all(link_df$gene %in% rownames(seur[[rna_assay]])))

#   Sample linked peaks from top, bottom, and middle ranges of FDR by cell type
link_df_list = list()
for (this_percentile in c(0, 0.5)) {
    link_df_list[[as.character(this_percentile)]] = link_df |>
        group_by(cell_type) |>
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
    group_by(cell_type) |>
    filter(n() > num_links * 4) |>
    arrange(desc(FDR)) |>
    slice_head(n = num_links) |>
    ungroup() |>
    mutate(percentile = 1)
link_sub_df = bind_rows(link_df_list)

count_df_list = list()
for (this_cell_type in unique(link_sub_df$cell_type)) {
    this_link_sub_df = link_sub_df |>
        filter(cell_type == this_cell_type)

    for (this_donor in unique(seur@meta.data$donor)) {
        subset_vec = (seur@meta.data$donor == this_donor) &
            (seur@meta.data$orig.ident == this_cell_type)
        
        count_df_list[[length(count_df_list) + 1]] = tibble(
            cell_type = this_link_sub_df$cell_type,
            donor = this_donor,
            percentile = this_link_sub_df$percentile,
            peak = this_link_sub_df$peak,
            gene = this_link_sub_df$gene,
            peak_value = rowMeans(
                GetAssayData(
                    seur, assay = atac_assay, layer = "data"
                )[this_link_sub_df$peak, subset_vec, drop = FALSE]
            ),
            gene_value = rowMeans(
                GetAssayData(
                    seur, assay = rna_assay, layer = "data"
                )[this_link_sub_df$gene, subset_vec, drop = FALSE]
            )
        )
    }
}
count_df = bind_rows(count_df_list) |>
    filter(!is.na(peak_value), !is.na(gene_value)) |>
    group_by(peak, gene, cell_type) |>
    #   Z-scoring across donors helps show the influence of donor independent of
    #   the magnitiudes of peak or gene expression
    mutate(
        peak_value = (peak_value - mean(peak_value)) / sd(peak_value),
        gene_value = (gene_value - mean(gene_value)) / sd(gene_value)
    ) |>
    ungroup()

cell_types <- unique(count_df$cell_type)
percentiles <- unique(count_df$percentile)

for (color_by in c("donor", "link")) {
    # Generate one plot per cell_type x percentile
    plots <- lapply(cell_types, function(ct) {
        lapply(percentiles, function(p) {
            df_sub <- count_df |>
                filter(cell_type == ct, percentile == p)
            
            if (color_by == 'link') {
                p = df_sub |>
                    mutate(link = paste(peak, gene, sep = "_")) |>
                    ggplot(aes(x = peak_value, y = gene_value, color = link)) +
                        scale_color_manual(values = link_colors)
            } else {
                p = ggplot(
                        df_sub,
                        aes(x = peak_value, y = gene_value, color = donor)
                    ) +
                    scale_color_manual(values = donor_colors)
            }

            p = p +
                geom_point(size = 0.5) +
                labs(
                    title = sprintf("%s | percentile=%.1f", ct, p),
                    x = "Peak", y = "Gene"
                ) +
                theme_bw(base_size = 15) +
                theme(legend.position = "none")
        })
    })

    # Flatten to a single list, row-major (cell_type changes slowly)
    plot_list <- unlist(plots, recursive = FALSE)

    # Extract shared legend
    if (color_by == 'link') {
        p = count_df |>
            mutate(link = paste(peak, gene, sep = "_")) |>
            ggplot(aes(x = peak_value, y = gene_value, color = link)) +
                scale_color_manual(values = link_colors)
    } else {
        p = ggplot(
                count_df,
                aes(x = peak_value, y = gene_value, color = donor)
            ) +
            scale_color_manual(values = donor_colors)
    }
    legend <- get_legend(
        p +
            geom_point() +
            theme_bw(base_size = 10) +
            guides(color = guide_legend(override.aes = list(size = 2)))
    )

    grid <- plot_grid(plotlist = plot_list, ncol = length(percentiles))
    final <- plot_grid(grid, legend, rel_widths = c(1, 0.15))

    pdf(
        file.path(plot_dir, sprintf("%s_%s_scatter.pdf", dataset, color_by)),
        width = 10, height = 3 * length(cell_types)
    )
    print(final)
    dev.off()
}


p = count_df |>
    group_by(cell_type, percentile, donor) |>
    #   Since we Z-scored across donors, a donor's typical distance from the
    #   origin can be an indirect measure of its contribution to the correlation
    summarize(
        contribution_score = mean(peak_value ** 2 + gene_value ** 2)
    ) |>
    ggplot(aes(x = donor, y = contribution_score, fill = donor)) +
        geom_bar(stat = "identity") +
        facet_grid(cell_type ~ percentile) +
        scale_fill_manual(values = donor_colors) +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
pdf(
    file.path(plot_dir, sprintf("%s_contribution_score_barplots.pdf", dataset)),
    width = 10, height = 1.5 * length(cell_types)
)
print(p)
dev.off()

p = count_df |>
    group_by(cell_type, donor) |>
    summarize(
        contribution_score = mean(peak_value ** 2 + gene_value ** 2)
    ) |>
    ggplot(aes(x = donor, y = contribution_score, color = donor)) +
        geom_boxplot() +
        scale_color_manual(values = donor_colors) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(x = "Donor", y = "Contribution Score") +
        guides(color = "none")
pdf(file.path(plot_dir, sprintf("%s_contribution_score_boxplots.pdf", dataset)))
print(p)
dev.off()

#   Compute a metric measuring the lopsidedness of donor contributions to
#   peak-gene correlations (higher is worse)
overall_score = count_df |>
    group_by(cell_type, donor) |>
    summarize(
        contribution_score = mean(peak_value ** 2 + gene_value ** 2)
    ) |>
    pull(contribution_score) |>
    var()
message(sprintf("Overall donor-specificity score: %.2f", overall_score))

session_info()
    