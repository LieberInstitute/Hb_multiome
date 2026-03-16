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

# Prepare matrices for each cell type
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

# Combine all matrices horizontally
combined_mat <- cbind(gene_mat_ct1, gene_mat_ct2, peak_mat_ct1, peak_mat_ct2)

anno_df <- count_df |>
    select(row_id, target_cell_type) |>
    distinct()

# Create color functions
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
    c("#440154FF", "#FDE725FF")
)

col_fun_peak <- colorRamp2(
    c(
        quantile(all_peak_values[all_peak_values < 0], 1 - color_quantile),
        quantile(all_peak_values[all_peak_values > 0], color_quantile)
    ), 
    c("#000004FF", "#FCFDBFFF")
)

cell_type_colors <- c(
    "LHb.2.7" = "#305252", "MHb.2" = "#F9B9B7", "both" = "#F5D491"
)

donors <- unique(count_df$donor)
donor_colors <- setNames(
    c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
      "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"),
    donors
)

# Create column split and data type vectors
col_split <- factor(
    c(
        rep(paste0("Gene: ", cell_type1), ncol(gene_mat_ct1)),
        rep(paste0("Gene: ", cell_type2), ncol(gene_mat_ct2)),
        rep(paste0("Peak: ", cell_type1), ncol(peak_mat_ct1)),
        rep(paste0("Peak: ", cell_type2), ncol(peak_mat_ct2))
    ),
    levels = c(
        paste0("Gene: ", cell_type1),
        paste0("Gene: ", cell_type2),
        paste0("Peak: ", cell_type1),
        paste0("Peak: ", cell_type2)
    )
)

data_type_vec <- c(
    rep("Gene", ncol(gene_mat_ct1) + ncol(gene_mat_ct2)),
    rep("Peak", ncol(peak_mat_ct1) + ncol(peak_mat_ct2))
)

cell_type_vec <- c(
    rep(cell_type1, ncol(gene_mat_ct1)),
    rep(cell_type2, ncol(gene_mat_ct2)),
    rep(cell_type1, ncol(peak_mat_ct1)),
    rep(cell_type2, ncol(peak_mat_ct2))
)

donor_vec <- c(
    colnames(gene_mat_ct1),
    colnames(gene_mat_ct2),
    colnames(peak_mat_ct1),
    colnames(peak_mat_ct2)
)

# Create column annotation
col_ha <- HeatmapAnnotation(
    `Cell Type` = cell_type_vec,
    Donor = donor_vec,
    col = list(`Cell Type` = cell_type_colors, Donor = donor_colors),
    show_legend = c(`Cell Type` = TRUE, Donor = FALSE)
)

# Split rows by target_cell_type
rows_ct1 <- which(anno_df$target_cell_type == cell_type1)
rows_ct2 <- which(anno_df$target_cell_type == cell_type2)
rows_both <- which(anno_df$target_cell_type == "both")

# Get row orderings by clustering dummy heatmaps on appropriate gene matrices
dummy_ht_ct1 <- Heatmap(gene_mat_ct1[rows_ct1, ], cluster_columns = FALSE)
row_order_ct1 <- row_order(dummy_ht_ct1)

dummy_ht_ct2 <- Heatmap(gene_mat_ct2[rows_ct2, ], cluster_columns = FALSE)
row_order_ct2 <- row_order(dummy_ht_ct2)

dummy_ht_both <- Heatmap(gene_mat_ct1[rows_both, ], cluster_columns = FALSE)
row_order_both <- row_order(dummy_ht_both)

# Combine row orders into a single vector
full_row_order <- c(rows_ct1[row_order_ct1], rows_ct2[row_order_ct2], rows_both[row_order_both])

# Create custom color function that handles both gene and peak values
col_fun_combined <- function(x, data_type) {
    if (data_type == "Gene") {
        col_fun_gene(x)
    } else {
        col_fun_peak(x)
    }
}

# Map colors to the combined matrix
color_mat <- matrix(NA, nrow = nrow(combined_mat), ncol = ncol(combined_mat))
for (i in 1:ncol(combined_mat)) {
    color_mat[, i] <- col_fun_combined(combined_mat[, i], data_type_vec[i])
}

# Create row annotation
row_ha <- rowAnnotation(
    `Target Cell Type` = anno_df$target_cell_type,
    col = list(`Target Cell Type` = cell_type_colors)
)

# Create the main heatmap
ht <- Heatmap(
    combined_mat,
    name = "Z-score",
    col = function(x) {
        col_fun_gene(x)
    },
    column_split = col_split,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_order = full_row_order,
    row_split = anno_df$target_cell_type,
    cluster_row_slices = FALSE,
    show_row_names = FALSE,
    left_annotation = row_ha,
    top_annotation = col_ha,
    border = TRUE,
    show_heatmap_legend = FALSE,  # Turn off automatic legend
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height,
                 gp = gpar(fill = color_mat[i, j], col = NA))
    }
)

# Create manual legends
gene_legend <- Legend(
    col_fun = col_fun_gene,
    title = "Gene Z-score",
    direction = "vertical"
)

peak_legend <- Legend(
    col_fun = col_fun_peak,
    title = "Peak Z-score",
    direction = "vertical"
)

pdf(plot_path, width = 10, height = 7)
draw(ht, annotation_legend_list = list(gene_legend, peak_legend))
dev.off()

session_info()
