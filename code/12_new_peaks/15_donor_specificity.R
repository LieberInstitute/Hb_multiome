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
    'processed-data', '12_new_peaks', '17_filter_links',
    sprintf('%s_filtered_data.parquet', dataset)
)
out_path = here(
    'processed-data', '12_new_peaks', '15_donor_specificity',
    sprintf('%s_filtered_unbiased.parquet', dataset)
)
plot_dir = here("plots", "12_new_peaks", "15_donor_specificity")
atac_assay = "ATAC"
rna_assay = "RNA"
num_expected_donors = 10
num_links_scatter = 20

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)
set.seed(0)

seur = qs_read(seur_path)
seur@meta.data$donor = str_extract(colnames(seur), 'S[0-9]{2}.*')
stopifnot(length(unique(seur@meta.data$donor)) == num_expected_donors)

donor_colors = palette36.colors(num_expected_donors)
names(donor_colors) = unique(seur@meta.data$donor)
link_colors = unname(palette36.colors(num_links_scatter))

link_df = read_parquet_duckdb(link_path) |>
    collect()

stopifnot(all(link_df$peak %in% rownames(seur[[atac_assay]])))
stopifnot(all(link_df$gene %in% rownames(seur[[rna_assay]])))

count_df_list = list()
for (this_cell_type in unique(link_df$cell_type)) {
    this_link_df = link_df |>
        filter(cell_type == this_cell_type)

    for (this_donor in unique(seur@meta.data$donor)) {
        subset_vec = (seur@meta.data$donor == this_donor) &
            (seur@meta.data$orig.ident == this_cell_type)
        
        count_df_list[[length(count_df_list) + 1]] = tibble(
            cell_type = this_link_df$cell_type,
            donor = this_donor,
            peak = this_link_df$peak,
            gene = this_link_df$gene,
            peak_value = rowMeans(
                GetAssayData(
                    seur, assay = atac_assay, layer = "data"
                )[this_link_df$peak, subset_vec, drop = FALSE]
            ),
            gene_value = rowMeans(
                GetAssayData(
                    seur, assay = rna_assay, layer = "data"
                )[this_link_df$gene, subset_vec, drop = FALSE]
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
    ungroup() |>
    group_by(cell_type, peak, gene) |>
    mutate(full_cor = cor(peak_value, gene_value)) |>
    #   For each peak-gene-cell_type, compute full-donor correlation, then
    #   leave-one-out (LOO) correlations. The donor-bias score is the drop in
    #   correlation when that donor is excluded.
    group_by(cell_type, peak, gene) |>
    mutate(full_cor = cor(peak_value, gene_value)) |>
    ungroup() |>
    #   For each row, compute correlation leaving out that donor
    group_by(cell_type, peak, gene) |>
    mutate(
        loo_cor = {
            x = peak_value
            y = gene_value
            n = length(x)
            vapply(seq_len(n), function(i) {
                if (n - 1 < 3) return(NA_real_)
                cor(x[-i], y[-i])
            }, numeric(1))
        },
        donor_bias = full_cor - loo_cor,
        #   The ratio of magnitudes between the larger and smaller correlations.
        #   The idea is we'll drop things with a value in this metric >= 1
        donor_score = case_when(
            (full_cor * loo_cor <= 0) | is.na(loo_cor) ~ 1,
            abs(full_cor) > abs(loo_cor) ~ full_cor / loo_cor - 1,
            TRUE ~ loo_cor / full_cor - 1
        )
    ) |>
    ungroup()

summary_df = count_df |>
    group_by(peak, gene, cell_type) |>
    summarize(donor_score = max(donor_score)) |>
    ungroup() |>
    left_join(link_df, by = c("peak", "gene", "cell_type"))

#   Check the distribution of donor scores
p = summary_df |>
    filter(donor_score != 1) |>
    ggplot(aes(x = donor_score)) +
    geom_histogram(
        bins = 60, fill = "steelblue", color = "white", linewidth = 0.2
    ) +
    scale_x_log10() +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    labs(x = "Donor Score (log10)", y = "Count") +
    theme_bw(base_size = 15)
pdf(
    file.path(plot_dir, sprintf("%s_donor_score_histogram.pdf", dataset)),
    height = 5
)
print(p)
dev.off()

message(
    sprintf(
        "Dropping %.1f%% of links with donor_score >= 1",
        mean(summary_df$donor_score >= 1) * 100
    )
)

#   Export the filtered set of links, dropping donor-biased ones
summary_df |>
    filter(donor_score < 1) |>
    select(all_of(colnames(link_df))) |>
    group_by(peak, gene) |>
    mutate(is_shared = n() > 1) |>
    ungroup() |>
    compute_parquet(out_path)

sample_df = count_df |>
    left_join(link_df, by = c("peak", "gene", "cell_type")) |>
    filter(score > 0) |>
    distinct(cell_type, peak, gene) |>
    group_by(cell_type) |>
    slice_sample(n = num_links_scatter) |>
    ungroup()

p = count_df |>
    inner_join(sample_df, by = c("cell_type", "peak", "gene")) |>
    ggplot(aes(x = peak_value, y = gene_value, color = donor)) +
        scale_color_manual(values = donor_colors) +
        geom_point(size = 0.5) +
        facet_wrap(~cell_type) +
        labs(x = "Peak", y = "Gene") +
        theme_bw(base_size = 20) +
        guides(color = guide_legend(override.aes = list(size = 4)))
pdf(
    file.path(plot_dir, sprintf("%s_scatter.pdf", dataset))
)
print(p)
dev.off()

bias_summary_df = count_df |>
    group_by(cell_type, donor) |>
    summarize(
        mean_donor_bias = mean(donor_bias, na.rm = TRUE), .groups = "drop"
    )

p_bias_bar = bias_summary_df |>
    ggplot(aes(x = donor, y = mean_donor_bias, fill = donor)) +
        geom_bar(stat = "identity") +
        facet_wrap(~cell_type) +
        scale_fill_manual(values = donor_colors) +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(x = "Donor", y = "Mean LOO Donor Bias")
pdf(file.path(plot_dir, sprintf("%s_loo_bias_barplots.pdf", dataset)))
print(p_bias_bar)
dev.off()

p_bias_box = bias_summary_df |>
    ggplot(aes(x = donor, y = mean_donor_bias, color = donor)) +
        geom_boxplot() +
        scale_color_manual(values = donor_colors) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        labs(x = "Donor", y = "Mean LOO Donor Bias") +
        guides(color = "none")
pdf(file.path(plot_dir, sprintf("%s_loo_bias_boxplots.pdf", dataset)))
print(p_bias_box)
dev.off()

#   Average (across all peak-gene pairs) of the largest absolute LOO bias
#   seen for any single donor. Interpretable as: "the most impactful donor
#   typically shifts the correlation by this much."
max_loo_score = count_df |>
    group_by(cell_type, peak, gene) |>
    summarize(
        max_abs_bias = max(abs(donor_bias), na.rm = TRUE), .groups = "drop"
    ) |>
    pull(max_abs_bias) |>
    mean(na.rm = TRUE)
message(
    sprintf(
        "Mean max-donor LOO bias (typical influence of most impactful donor): %.3f",
        max_loo_score
    )
)

session_info()
