library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(ComplexHeatmap)
library(qs2)
library(Seurat)

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)

prep_path = here(
    "processed-data", "13_tripod_trios", "02_tripod_preprocess",
    "preprocessed_objects_%s.qs2"
)

plot_dir = here("plots", "13_tripod_trios", "13_summary_heatmap")

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

#-------------------------------------------------------------------------------
#   Load and filter trio data
#-------------------------------------------------------------------------------

# Load top 20 trios per cell type by adjusted p-value (stringency level 1)
trio_df = read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    filter(is_intersect, stringency_level == 1) |>
    collect() |>
    group_by(cell_type) |>
    arrange(adj) |>
    slice(1:20) |>
    ungroup() |>
    select(peak, gene, TF, cell_type, adj)

#-------------------------------------------------------------------------------
#   Load metacell data and build feature matrices
#-------------------------------------------------------------------------------

all_cts <- unique(trio_df$cell_type)
mats_list <- setNames(vector("list", length(all_cts)), all_cts)

for (ct in all_cts) {
    mats_list[[ct]] <- qs_read(sprintf(prep_path, ct))$metacell_seur
}

# Determine cell type order dynamically from the data to support subsetting
# Order will be: cell types as they appear in the data, in a consistent order
cell_type_order <- as.character(unique(factor(trio_df$cell_type, levels = unique(trio_df$cell_type))))

# Build data matrices for TF, peak, and gene across all trios and metacells
tf_mat <- list()
peak_mat <- list()
gene_mat <- list()

for (i in seq_len(nrow(trio_df))) {
    row <- trio_df[i, ]
    tf_name <- row$TF
    peak_name <- row$peak
    gene_name <- row$gene
    
    tf_vals <- c()
    peak_vals <- c()
    gene_vals <- c()
    col_names <- c()
    
    # Iterate through cell types in consistent order
    for (ct in cell_type_order) {
        mats <- mats_list[[ct]]
        n_metacells <- nrow(mats$rna)
        
        # Get values for this feature; assume 0 if feature not in matrix
        tf_row <- tryCatch(
            mats$rna[, tf_name, drop = TRUE],
            error = function(e) rep(0, n_metacells)
        )
        if (length(tf_row) == 0) tf_row <- rep(0, n_metacells)
        
        gene_row <- tryCatch(
            mats$rna[, gene_name, drop = TRUE],
            error = function(e) rep(0, n_metacells)
        )
        if (length(gene_row) == 0) gene_row <- rep(0, n_metacells)
        
        peak_row <- tryCatch(
            mats$peak[, peak_name, drop = TRUE],
            error = function(e) rep(0, n_metacells)
        )
        if (length(peak_row) == 0) peak_row <- rep(0, n_metacells)
        
        tf_vals <- c(tf_vals, tf_row)
        gene_vals <- c(gene_vals, gene_row)
        peak_vals <- c(peak_vals, peak_row)
        
        # Create column names: ct_metacell_0, ct_metacell_1, etc
        metacell_names <- paste0(ct, "_", seq_len(n_metacells) - 1)
        col_names <- c(col_names, metacell_names)
    }
    
    tf_mat[[i]] <- tf_vals
    peak_mat[[i]] <- peak_vals
    gene_mat[[i]] <- gene_vals
    
    if (i == 1) {
        stored_col_names <- col_names
    }
}

# Convert to matrices and Z-score each row independently
tf_mat <- do.call(rbind, tf_mat)
peak_mat <- do.call(rbind, peak_mat)
gene_mat <- do.call(rbind, gene_mat)

colnames(tf_mat) <- stored_col_names
colnames(peak_mat) <- stored_col_names
colnames(gene_mat) <- stored_col_names

z_score_row <- function(row) {
    m <- mean(row)
    s <- sd(row)
    if (s == 0) return(row - m)
    return((row - m) / s)
}

tf_mat_z <- t(apply(tf_mat, 1, z_score_row))
peak_mat_z <- t(apply(peak_mat, 1, z_score_row))
gene_mat_z <- t(apply(gene_mat, 1, z_score_row))

# Combine matrices: TF, Peak, Gene (rows will be ordered by cell type then significance)
combined_mat <- rbind(tf_mat_z, peak_mat_z, gene_mat_z)
colnames(combined_mat) <- stored_col_names

#-------------------------------------------------------------------------------
#   Create ComplexHeatmap with annotations
#-------------------------------------------------------------------------------

# Row annotations
row_groups <- c(
    rep("TF", nrow(tf_mat_z)),
    rep("Peak", nrow(peak_mat_z)),
    rep("Gene", nrow(gene_mat_z))
)

row_names_tf <- paste0("TF: ", trio_df$TF)
row_names_peak <- paste0("Peak: ", substr(trio_df$peak, 1, 20))
row_names_gene <- paste0("Gene: ", trio_df$gene)
row_names_combined <- c(row_names_tf, row_names_peak, row_names_gene)

rownames(combined_mat) <- row_names_combined

row_annot <- rowAnnotation(
    feature_type = row_groups,
    col = list(feature_type = c("TF" = "lightblue", "Peak" = "lightgreen", "Gene" = "lightyellow")),
    show_legend = TRUE,
    gp = gpar(fontsize = 10)
)

# Column annotations
col_cell_types <- sub("_[0-9]+$", "", colnames(combined_mat))

col_annot <- HeatmapAnnotation(
    cell_type = col_cell_types,
    col = list(cell_type = setNames(
        scales::hue_pal()(length(unique(col_cell_types))),
        unique(col_cell_types)
    )),
    show_legend = TRUE,
    gp = gpar(fontsize = 8),
    simple_anno_size = unit(0.5, "cm")
)

# Create heatmap
hm <- Heatmap(
    combined_mat,
    name = "Z-score",
    col = circlize::colorRamp2(
        breaks = c(-2, 0, 2),
        colors = c("blue", "white", "red")
    ),
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    cluster_rows = FALSE,  # Rows ordered by cell type then significance
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    left_annotation = row_annot,
    top_annotation = col_annot,
    row_split = row_groups,
    column_split = col_cell_types,
    width = unit(22, "cm"),
    height = unit(25, "cm"),
    use_raster = TRUE,
    raster_quality = 2,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(
        title = "Z-score",
        title_gp = gpar(fontsize = 11, fontface = "bold"),
        at = c(-2, 0, 2),
        direction = "vertical",
        legend_width = unit(2, "cm"),
        legend_height = unit(4, "cm")
    )
)

# Save heatmap
pdf(file.path(plot_dir, "summary_heatmap.pdf"), width = 14, height = 16)
draw(
    hm,
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "bottom",
    padding = unit(c(2, 8, 2, 2), "mm")
)
dev.off()

session_info()
