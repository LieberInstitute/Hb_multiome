#   How much are overall peak-gene correlations measured by linked peaks
#   influenced by a handful of donors? Produce various visualizations and
#   metrics to quantify this bias (for metacell-derived links in LHb.2.7,
#   as opposed to donor-pseudobulked links from all cell types)

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
    "processed-data", "12_new_peaks", "14_metacell_link_peaks",
    "LHb.2.7_LHb.2.7.csv.gz"
)
plot_dir = here("plots", "12_new_peaks", "18_donor_specificity_metacell")
atac_assay = "ATAC_macs2_pseudo"
rna_assay = "RNA"
this_cell_type = "LHb.2.7"
num_expected_donors = 10
num_links = 20
cor_thres = 0.2
FDR_thres = 0.2

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

seur = readRDS(seur_path)
seur@meta.data$donor = str_extract(colnames(seur), 'S[0-9]{2}.*')
stopifnot(length(unique(seur@meta.data$donor)) == num_expected_donors)

donor_colors = palette36.colors(num_expected_donors)
names(donor_colors) = unique(seur@meta.data$donor)

link_df = read_csv_duckdb(link_path, prudence = "lavish") |>
    filter(
        FDR < FDR_thres, score > cor_thres,
        gene %in% rownames(seur[[rna_assay]]),
        peak %in% rownames(seur[[atac_assay]])
    ) |>
    select(peak, gene, other_cell_type, score, FDR) |>
    rename(cell_type = other_cell_type) |>
    collect()

#   Sample linked peaks from top, bottom, and middle ranges of FDR
link_df_list = list()
for (this_percentile in c(0, 0.5)) {
    link_df_list[[as.character(this_percentile)]] = link_df |>
        filter(
            FDR > quantile(FDR, this_percentile),
            n() > num_links * 4
        ) |>
        arrange(FDR) |>
        slice_head(n = num_links) |>
        mutate(percentile = this_percentile)
}
link_df_list[['1']] = link_df |>
    filter(n() > num_links * 4) |>
    arrange(desc(FDR)) |>
    slice_head(n = num_links) |>
    mutate(percentile = 1)
link_sub_df = bind_rows(link_df_list)

count_df_list = list()
for (this_donor in unique(seur@meta.data$donor)) {
    subset_vec = (seur@meta.data$donor == this_donor) &
        (seur@meta.data$orig.ident == this_cell_type)
    
    count_df_list[[length(count_df_list) + 1]] = tibble(
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
    filter(!is.na(peak_value), !is.na(gene_value)) |>
    group_by(peak, gene) |>
    #   Z-scoring across donors helps show the influence of donor independent of
    #   the magnitiudes of peak or gene expression
    mutate(
        peak_value = (peak_value - mean(peak_value)) / sd(peak_value),
        gene_value = (gene_value - mean(gene_value)) / sd(gene_value)
    ) |>
    ungroup()

p = ggplot(count_df, aes(x = peak_value, y = gene_value, color = donor)) +
    geom_point(size = 0.5) +
    scale_color_manual(values = donor_colors) +
    facet_wrap(~ percentile, nrow = 1) +
    labs(x = "Peak", y = "Gene") +
    theme_bw(base_size = 20)
pdf(
    file.path(plot_dir, "specificity_scatter.pdf"),
    width = 12, height = 4
)
print(p)
dev.off()

p = count_df |>
    group_by(percentile, donor) |>
    #   Since we Z-scored across donors, a donor's typical distance from the
    #   origin can be an indirect measure of its contribution to the correlation
    summarize(
        contribution_score = mean(peak_value ** 2 + gene_value ** 2)
    ) |>
    ggplot(aes(x = donor, y = contribution_score, fill = donor)) +
        geom_bar(stat = "identity") +
        facet_wrap(~ percentile, nrow = 3) +
        scale_fill_manual(values = donor_colors) +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
pdf(
    file.path(plot_dir, "contribution_score_barplots.pdf")
)
print(p)
dev.off()

#   Compute a metric measuring the lopsidedness of donor contributions to
#   peak-gene correlations (higher is worse)
overall_score = count_df |>
    group_by(donor) |>
    summarize(
        contribution_score = mean(peak_value ** 2 + gene_value ** 2)
    ) |>
    pull(contribution_score) |>
    var()
message(sprintf("Overall donor-specificity score: %.2f", overall_score))

session_info()
    