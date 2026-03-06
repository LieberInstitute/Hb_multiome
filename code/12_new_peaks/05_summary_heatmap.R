library(sessioninfo)
library(Seurat)
library(Signac)
library(GenomicRanges) 
library(tidyverse)
library(here)
library(Matrix)
library(ComplexHeatmap)
library(circlize)

seur_path = here(
    "processed-data", "06_peak_calling", "12_pseudobulk_MACS2",
    "Mid_pseudobulk.spearman.5e5.rds"
)
result_paths = here(
    'processed-data', '12_new_peaks', '01_link_peaks', '%s_%s.csv.gz'
)
plot_path = here(
    "plots", "12_new_peaks", "05_summary_heatmap", "heatmap.pdf"
)
atac_assay = "ATAC_macs2_pseudo"
rna_assay = "RNA"
cell_type1 = "MHb.2"
cell_type2 = "LHb.2.7"
cor_thres = 0.3
FDR_thres = 0.1
num_expected_donors = 10
color_quantile = 0.99

dir.create(dirname(plot_path), showWarnings = FALSE)

seur = readRDS(seur_path)
seur@meta.data$donor = str_extract(colnames(seur), 'S[0-9]{2}.*')

#   We want a heatmap without NA values, so check all donors are present in the
#   relevant cell types
stopifnot(sum(seur@meta.data$orig.ident %in% c(cell_type1, cell_type2)) == 2 * num_expected_donors)

link_df = rbind(
        read_csv(
            sprintf(result_paths, cell_type1, cell_type1),
            show_col_types = FALSE
        ),
        read_csv(
            sprintf(result_paths, cell_type2, cell_type2),
            show_col_types = FALSE
        )
    ) |>
    filter(abs(score) > cor_thres, FDR < FDR_thres) |>
    select(peak, gene, target_cell_type, score) |>
    group_by(peak, gene) |>
    summarize(
        target_cell_type = ifelse(n() == 1, target_cell_type, "both"),
        score = first(score) # this is ok since scores will be identical
    ) |>
    ungroup()

stopifnot(all(link_df$peak %in% rownames(seur[[atac_assay]])))
stopifnot(all(link_df$gene %in% rownames(seur[[rna_assay]])))

count_df_list = list()
for (cell_type in unique(link_df$target_cell_type)) {
    # Subset links for this cell type
    ct_links <- link_df |>
        filter(target_cell_type == cell_type)
  
    for (donor in unique(seur@meta.data$donor)) {
        #   Determine how to subset by cell type and donor
        if (cell_type == "both") {
            subset_vec = (
                (seur@meta.data$orig.ident %in% c(cell_type1, cell_type2)) &
                (seur@meta.data$donor == donor)
            )
        } else {
            subset_vec = (
                (seur@meta.data$orig.ident == cell_type) &
                (seur@meta.data$donor == donor)
            )
        }
      
        count_df_list[[length(count_df_list) + 1]] <- tibble(
            cell_type = cell_type,
            donor = donor,
            peak = ct_links$peak,
            gene = ct_links$gene,
            score = ct_links$score,
            peak_value = rowMeans(
                GetAssayData(
                    seur, assay = atac_assay, layer = "data"
                )[ct_links$peak, subset_vec, drop = FALSE]
            ),
            gene_value = rowMeans(
                GetAssayData(
                    seur, assay = rna_assay, layer = "data"
                )[ct_links$gene, subset_vec, drop = FALSE]
            )
        )
    }
}
count_df = bind_rows(count_df_list) |>
    mutate(row_id = paste(peak, gene, sep = "|")) |>
    group_by(row_id) |>
    mutate(
        peak_value = (peak_value - mean(peak_value)) / sd(peak_value),
        gene_value = (gene_value - mean(gene_value)) / sd(gene_value)
    ) |>
    ungroup()

# Prepare matrices for peak and gene values
peak_mat <- count_df |>
    select(row_id, donor, peak_value) |>
    pivot_wider(names_from = donor, values_from = peak_value) |>
    column_to_rownames("row_id") |>
    as.matrix()

gene_mat <- count_df |>
    select(row_id, donor, gene_value) |>
    pivot_wider(names_from = donor, values_from = gene_value) |>
    column_to_rownames("row_id") |>
    as.matrix()
stopifnot(identical(rownames(peak_mat), rownames(gene_mat)))

# Create row split factor with cell type grouped together
row_split_df <- count_df |>
    select(peak, gene, cell_type, score) |>
    distinct() |>
    mutate(
        cor_sign = ifelse(score > 0, "positive", "negative"),
        split_group = paste(cell_type, cor_sign, sep = "_")
    )
row_split_levels = row_split_df |>
    distinct(cell_type, cor_sign) |>
    arrange(cell_type, cor_sign) |>
    mutate(split_group = paste(cell_type, cor_sign, sep = "_")) |>
    pull(split_group)
row_split = factor(row_split_df$split_group, levels = row_split_levels)

# Create annotation data frame with same ordering
anno_df <- count_df |>
    select(peak, gene, cell_type, score) |>
    distinct() |>
    mutate(cor_sign = ifelse(score > 0, "positive", "negative"))

# Create color function
col_fun_gene <- colorRamp2(
    c(
        quantile(gene_mat[gene_mat < 0], 1 - color_quantile),
        quantile(gene_mat[gene_mat > 0], color_quantile)
    ),
    c("#440154FF", "#FDE725FF")  # viridis colors
)

col_fun_peak <- colorRamp2(
    c(
        quantile(peak_mat[peak_mat < 0], 1 - color_quantile),
        quantile(peak_mat[peak_mat > 0], color_quantile)
    ), 
    c("#000004FF", "#FCFDBFFF")  # magma colors
)

# Define colors for annotations
cell_type_colors <- c(
    "LHb.2.7" = "#305252", "MHb.2" = "#F9B9B7", "both" = "#F5D491"
)
cor_sign_colors <- c("positive" = "#6D9DC5", "negative" = "#E3170A")

# Create row annotation
row_ha <- rowAnnotation(
    `Cell Type` = factor(anno_df$cell_type, levels = names(cell_type_colors)),
    `Correlation` = factor(anno_df$cor_sign, levels = names(cor_sign_colors)),
    col = list(
        `Cell Type` = cell_type_colors,
        `Correlation` = cor_sign_colors
    ),
    show_legend = TRUE
)

#   RNA heatmap
ht_gene <- Heatmap(
    gene_mat,
    name = "Gene Z-score",
    col = col_fun_gene,
    column_title = "Gene Values",
    show_row_names = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    row_split = row_split,
    cluster_row_slices = FALSE,
    row_title = NULL,  # Remove row split labels
    left_annotation = row_ha,
    border = TRUE
)

# Get the row order from the gene heatmap
gene_row_order <- row_order(ht_gene)

#   ATAC heatmap matching row order
ht_peak <- Heatmap(
    peak_mat,
    name = "Peak Z-score",
    col = col_fun_peak,
    column_title = "Peak Values",
    show_row_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    row_split = row_split,
    row_title = NULL,  # Remove row split labels
    row_order = unlist(gene_row_order),
    border = TRUE
)

# Combine heatmaps
ht_list <- ht_gene + ht_peak
pdf(plot_path)
draw(ht_list)
dev.off()
