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
    #   Note the lack of abs() is intentional here; we only want to visualize
    #   positive correlations in the heatmap
    filter(score > cor_thres, FDR < FDR_thres) |>
    select(peak, gene, target_cell_type) |>
    group_by(peak, gene) |>
    summarize(
        target_cell_type = ifelse(n() == 1, unique(target_cell_type), "both")
    ) |>
    ungroup()

stopifnot(all(link_df$peak %in% rownames(seur[[atac_assay]])))
stopifnot(all(link_df$gene %in% rownames(seur[[rna_assay]])))

count_df_list = list()
for (cell_type in c(cell_type1, cell_type2)) {
    for (donor in unique(seur@meta.data$donor)) {
        #   Determine how to subset by cell type and donor
        subset_vec = (
            (seur@meta.data$orig.ident == cell_type) &
            (seur@meta.data$donor == donor)
        )
      
        count_df_list[[length(count_df_list) + 1]] <- tibble(
            cell_type = cell_type,
            target_cell_type = link_df$target_cell_type,
            donor = donor,
            peak = link_df$peak,
            gene = link_df$gene,
            peak_value = rowMeans(
                GetAssayData(
                    seur, assay = atac_assay, layer = "data"
                )[link_df$peak, subset_vec, drop = FALSE]
            ),
            gene_value = rowMeans(
                GetAssayData(
                    seur, assay = rna_assay, layer = "data"
                )[link_df$gene, subset_vec, drop = FALSE]
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
    ungroup() |>
    mutate(
        target_cell_type = factor(
            target_cell_type, levels = c(cell_type1, cell_type2, "both")
        )
    ) |>
    arrange(target_cell_type)

# Prepare matrices for each cell type separately
gene_mat_ct1 <- count_df |>
    filter(cell_type == cell_type1) |>
    select(row_id, donor, gene_value) |>
    pivot_wider(names_from = donor, values_from = gene_value) |>
    column_to_rownames("row_id") |>
    as.matrix()

gene_mat_ct2 <- count_df |>
    filter(cell_type == cell_type2) |>
    select(row_id, donor, gene_value) |>
    pivot_wider(names_from = donor, values_from = gene_value) |>
    column_to_rownames("row_id") |>
    as.matrix()

peak_mat_ct1 <- count_df |>
    filter(cell_type == cell_type1) |>
    select(row_id, donor, peak_value) |>
    pivot_wider(names_from = donor, values_from = peak_value) |>
    column_to_rownames("row_id") |>
    as.matrix()

peak_mat_ct2 <- count_df |>
    filter(cell_type == cell_type2) |>
    select(row_id, donor, peak_value) |>
    pivot_wider(names_from = donor, values_from = peak_value) |>
    column_to_rownames("row_id") |>
    as.matrix()

# Create row split factor based on target_cell_type
row_split_df <- count_df |>
    select(row_id, target_cell_type) |>
    distinct()
row_split <- row_split_df$target_cell_type

# Create annotation data frame
anno_df <- count_df |>
    select(row_id, target_cell_type) |>
    distinct()

# Create color function
all_gene_values <- c(
    count_df$gene_value[count_df$gene_value < 0],
    count_df$gene_value[count_df$gene_value > 0]
)
all_peak_values <- c(
    count_df$peak_value[count_df$peak_value < 0],
    count_df$peak_value[count_df$peak_value > 0]
)

col_fun_gene <- colorRamp2(
    c(
        quantile(all_gene_values[all_gene_values < 0], 1 - color_quantile),
        quantile(all_gene_values[all_gene_values > 0], color_quantile)
    ),
    c("#440154FF", "#FDE725FF")  # viridis colors
)

col_fun_peak <- colorRamp2(
    c(
        quantile(all_peak_values[all_peak_values < 0], 1 - color_quantile),
        quantile(all_peak_values[all_peak_values > 0], color_quantile)
    ), 
    c("#000004FF", "#FCFDBFFF")  # magma colors
)

# Define colors for annotations
cell_type_colors <- c(
    "LHb.2.7" = "#305252", "MHb.2" = "#F9B9B7", "both" = "#F5D491"
)

# Create donor colors (will not show legend)
donors <- unique(count_df$donor)
donor_colors <- setNames(
    rainbow(length(donors)),
    donors
)

# Create row annotation
row_ha <- rowAnnotation(
    `Target Cell Type` = anno_df$target_cell_type,
    col = list(
        `Target Cell Type` = cell_type_colors
    ),
    show_legend = TRUE
)

# Create column annotations for each heatmap
col_ha_gene_ct1 <- HeatmapAnnotation(
    `Cell Type` = rep(cell_type1, ncol(gene_mat_ct1)),
    Donor = colnames(gene_mat_ct1),
    col = list(
        `Cell Type` = cell_type_colors,
        Donor = donor_colors
    ),
    show_legend = c(`Cell Type` = TRUE, Donor = FALSE)
)

col_ha_gene_ct2 <- HeatmapAnnotation(
    `Cell Type` = rep(cell_type2, ncol(gene_mat_ct2)),
    Donor = colnames(gene_mat_ct2),
    col = list(
        `Cell Type` = cell_type_colors,
        Donor = donor_colors
    ),
    show_legend = c(`Cell Type` = FALSE, Donor = FALSE)
)

col_ha_peak_ct1 <- HeatmapAnnotation(
    `Cell Type` = rep(cell_type1, ncol(peak_mat_ct1)),
    Donor = colnames(peak_mat_ct1),
    col = list(
        `Cell Type` = cell_type_colors,
        Donor = donor_colors
    ),
    show_legend = c(`Cell Type` = FALSE, Donor = FALSE)
)

col_ha_peak_ct2 <- HeatmapAnnotation(
    `Cell Type` = rep(cell_type2, ncol(peak_mat_ct2)),
    Donor = colnames(peak_mat_ct2),
    col = list(
        `Cell Type` = cell_type_colors,
        Donor = donor_colors
    ),
    show_legend = c(`Cell Type` = FALSE, Donor = FALSE)
)

# Create heatmaps
ht_gene_ct1 <- Heatmap(
    gene_mat_ct1,
    name = "Gene Z-score",
    col = col_fun_gene,
    column_title = paste0("Gene: ", cell_type1),
    show_row_names = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    row_split = row_split,
    cluster_row_slices = FALSE,
    left_annotation = row_ha,
    top_annotation = col_ha_gene_ct1,
    border = TRUE
)

# Get row order from first heatmap
row_order_list <- row_order(ht_gene_ct1)

ht_gene_ct2 <- Heatmap(
    gene_mat_ct2,
    name = "Gene Z-score",
    col = col_fun_gene,
    column_title = paste0("Gene: ", cell_type2),
    show_row_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    row_split = row_split,
    row_order = unlist(row_order_list),
    top_annotation = col_ha_gene_ct2,
    border = TRUE,
    show_heatmap_legend = FALSE
)

ht_peak_ct1 <- Heatmap(
    peak_mat_ct1,
    name = "Peak Z-score",
    col = col_fun_peak,
    column_title = paste0("Peak: ", cell_type1),
    show_row_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    row_split = row_split,
    row_order = unlist(row_order_list),
    top_annotation = col_ha_peak_ct1,
    border = TRUE
)

ht_peak_ct2 <- Heatmap(
    peak_mat_ct2,
    name = "Peak Z-score",
    col = col_fun_peak,
    column_title = paste0("Peak: ", cell_type2),
    show_row_names = FALSE,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    row_split = row_split,
    row_order = unlist(row_order_list),
    top_annotation = col_ha_peak_ct2,
    border = TRUE,
    show_heatmap_legend = FALSE
)

# Combine heatmaps
ht_list <- ht_gene_ct1 + ht_gene_ct2 + ht_peak_ct1 + ht_peak_ct2
pdf(plot_path)
draw(ht_list)
dev.off()
