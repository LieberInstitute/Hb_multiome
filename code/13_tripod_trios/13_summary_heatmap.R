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
num_cores = if (is.na(num_cores)) 1L else num_cores
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

#-------------------------------------------------------------------------------
#   Load and filter trio data
#-------------------------------------------------------------------------------

# Top 20 trios per cell type by adjusted p-value (stringency level 1)
trio_df = read_parquet_duckdb(trio_path, prudence = "stingy") |>
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

cell_type_order <- unique(trio_df$cell_type)

mats_list <- setNames(vector("list", length(cell_type_order)), cell_type_order)
for (ct in cell_type_order) {
    mats_list[[ct]] <- qs_read(sprintf(prep_path, ct))$metacell_seur
}

# Helper: extract a feature vector from a matrix, returning zeros if absent
get_feature <- function(mat, name) {
    n <- nrow(mat)
    vals <- tryCatch(mat[, name, drop = TRUE], error = function(e) numeric(n))
    if (length(vals) == 0) numeric(n) else vals
}

# Raw feature matrices: rows = trios, cols = all metacells concatenated by CT
stored_col_names <- unlist(lapply(cell_type_order, function(ct) {
    paste0(ct, "_", seq_len(nrow(mats_list[[ct]]$rna)) - 1)
}))

tf_mat   <- do.call(rbind, lapply(seq_len(nrow(trio_df)), function(i) {
    unlist(lapply(cell_type_order, function(ct) get_feature(mats_list[[ct]]$rna, trio_df$TF[i])))
}))
peak_mat <- do.call(rbind, lapply(seq_len(nrow(trio_df)), function(i) {
    unlist(lapply(cell_type_order, function(ct) get_feature(mats_list[[ct]]$peak, trio_df$peak[i])))
}))
gene_mat <- do.call(rbind, lapply(seq_len(nrow(trio_df)), function(i) {
    unlist(lapply(cell_type_order, function(ct) get_feature(mats_list[[ct]]$rna, trio_df$gene[i])))
}))

colnames(tf_mat) <- colnames(peak_mat) <- colnames(gene_mat) <- stored_col_names

# Column index ranges for each cell type
ct_col_ranges <- local({
    start <- 1L
    lapply(setNames(cell_type_order, cell_type_order), function(ct) {
        n <- nrow(mats_list[[ct]]$rna)
        idx <- start:(start + n - 1L)
        start <<- start + n
        idx
    })
})

# Shared color palettes
ct_colors   <- setNames(scales::hue_pal()(length(cell_type_order)), cell_type_order)
feat_colors <- c("Gene" = "#4DAF4A", "Peak" = "#377EB8", "TF" = "#E41A1C")

#-------------------------------------------------------------------------------
#   Heatmap approach 1: rows = features (TF/Peak/Gene stacked), cols = metacells
#
#   Z-scoring: each row Z-scored globally across all metacells.
#   Row layout: three stacked blocks (TF, Peak, Gene), one row per trio per block.
#   Row ordering: fixed by cell type then significance (no clustering).
#   Column layout: one slice per cell type, metacells clustered within slice.
#-------------------------------------------------------------------------------

z_score_row <- function(row) {
    s <- sd(row)
    (row - mean(row)) / if (s == 0) 1 else s
}

tf_mat_z   <- t(apply(tf_mat,   1, z_score_row))
peak_mat_z <- t(apply(peak_mat, 1, z_score_row))
gene_mat_z <- t(apply(gene_mat, 1, z_score_row))

combined_mat <- rbind(tf_mat_z, peak_mat_z, gene_mat_z)
colnames(combined_mat) <- stored_col_names
rownames(combined_mat) <- c(
    paste0("TF: ",   trio_df$TF),
    paste0("Peak: ", substr(trio_df$peak, 1, 20)),
    paste0("Gene: ", trio_df$gene)
)

row_groups    <- rep(c("TF", "Peak", "Gene"), each = nrow(trio_df))
col_cell_types <- sub("_[0-9]+$", "", colnames(combined_mat))

hm1 <- Heatmap(
    combined_mat,
    name             = "Z-score",
    col              = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    cluster_rows     = FALSE,
    cluster_columns  = TRUE,
    cluster_column_slices = FALSE,
    show_row_dend    = FALSE,
    show_column_dend = FALSE,
    show_row_names   = FALSE,
    show_column_names = FALSE,
    row_split        = row_groups,
    column_split     = col_cell_types,
    left_annotation  = rowAnnotation(
        feature_type = row_groups,
        col = list(feature_type = c("TF" = "lightblue", "Peak" = "lightgreen", "Gene" = "lightyellow")),
        show_legend = TRUE,
        gp = gpar(fontsize = 10)
    ),
    top_annotation   = HeatmapAnnotation(
        cell_type = col_cell_types,
        col = list(cell_type = ct_colors),
        show_legend = TRUE,
        gp = gpar(fontsize = 8),
        simple_anno_size = unit(0.5, "cm")
    ),
    width            = unit(22, "cm"),
    height           = unit(25, "cm"),
    use_raster       = TRUE,
    raster_quality   = 2,
    column_title_gp  = gpar(fontsize = 10, fontface = "bold"),
    row_title_gp     = gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(
        title      = "Z-score",
        title_gp   = gpar(fontsize = 11, fontface = "bold"),
        at         = c(-2, 0, 2),
        direction  = "vertical",
        legend_height = unit(4, "cm")
    )
)

pdf(file.path(plot_dir, "summary_heatmap_approach1.pdf"), width = 14, height = 16)
draw(hm1, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "bottom", padding = unit(c(2, 8, 2, 2), "mm"))
dev.off()

#-------------------------------------------------------------------------------
#   Heatmap approach 2: rows = trios, cols = metacells grouped as Gene|Peak|TF
#
#   Z-scoring: each feature Z-scored within each cell type's metacell block,
#     applied uniformly regardless of which CT the trio was defined in.
#   Row layout: one row per trio (all three features encoded as column groups).
#   Row ordering: grouped by source CT; within each group, hierarchical
#     clustering on all three features from the source CT's metacell columns.
#   Column layout: for each CT in order, metacells appear three times
#     side-by-side — Gene | Peak | TF — in the same metacell order.
#-------------------------------------------------------------------------------

# Z-score each feature within each CT's metacell block independently
z_score_block <- function(mat, col_idx) {
    t(apply(mat[, col_idx, drop = FALSE], 1, z_score_row))
}

gene_mat_z2 <- gene_mat
peak_mat_z2 <- peak_mat
tf_mat_z2   <- tf_mat
for (ct in cell_type_order) {
    idx <- ct_col_ranges[[ct]]
    gene_mat_z2[, idx] <- z_score_block(gene_mat, idx)
    peak_mat_z2[, idx] <- z_score_block(peak_mat, idx)
    tf_mat_z2[,   idx] <- z_score_block(tf_mat,   idx)
}

# Build new column layout: Gene | Peak | TF metacells for each CT in sequence
col_blocks <- lapply(cell_type_order, function(ct) {
    idx <- ct_col_ranges[[ct]]
    mc  <- paste0(ct, "_", seq_along(idx) - 1)
    list(
        mat  = cbind(gene_mat_z2[, idx], peak_mat_z2[, idx], tf_mat_z2[, idx]),
        names = c(paste0(ct, "_gene_", mc), paste0(ct, "_peak_", mc), paste0(ct, "_tf_", mc)),
        ct   = rep(ct,  3 * length(idx)),
        feat = rep(c("Gene", "Peak", "TF"), each = length(idx))
    )
})

new_mat        <- do.call(cbind, lapply(col_blocks, `[[`, "mat"))
col_annot_ct   <- unlist(lapply(col_blocks, `[[`, "ct"))
col_annot_feat <- unlist(lapply(col_blocks, `[[`, "feat"))
colnames(new_mat) <- unlist(lapply(col_blocks, `[[`, "names"))

# Row ordering: group by source CT, cluster within each group using all three
# features from that CT's metacell columns
source_ct <- trio_df$cell_type
row_order  <- integer(0)
row_group_labels <- character(0)

for (ct in cell_type_order) {
    trio_idx <- which(source_ct == ct)
    if (length(trio_idx) == 0) next
    sub_mat  <- new_mat[trio_idx, col_annot_ct == ct, drop = FALSE]
    ordering <- if (length(trio_idx) == 1) 1L else hclust(dist(sub_mat), method = "complete")$order
    row_order        <- c(row_order, trio_idx[ordering])
    row_group_labels <- c(row_group_labels, rep(ct, length(trio_idx)))
}

new_mat_ordered <- new_mat[row_order, ]
row_source_ct   <- source_ct[row_order]
row_labels      <- paste0(trio_df$gene[row_order], " / ", trio_df$TF[row_order])

hm2 <- Heatmap(
    new_mat_ordered,
    name = "Z-score",
    col  = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
    cluster_columns       = FALSE,
    cluster_column_slices = FALSE,
    column_split          = factor(col_annot_ct, levels = cell_type_order),
    show_column_names     = FALSE,
    show_column_dend      = FALSE,
    column_title_gp       = gpar(fontsize = 7, fontface = "bold"),
    column_title_rot      = 90,
    column_gap            = unit(1.5, "mm"),
    cluster_rows          = FALSE,
    cluster_row_slices    = FALSE,
    row_split             = factor(row_group_labels, levels = cell_type_order),
    show_row_dend         = FALSE,
    show_row_names        = FALSE,
    row_title_gp          = gpar(fontsize = 7, fontface = "bold"),
    row_gap               = unit(1, "mm"),
    top_annotation = HeatmapAnnotation(
        `Cell type` = col_annot_ct,
        Feature     = col_annot_feat,
        col = list(`Cell type` = ct_colors, Feature = feat_colors),
        show_legend = TRUE,
        simple_anno_size = unit(0.4, "cm"),
        annotation_name_side = "left",
        annotation_name_gp   = gpar(fontsize = 8)
    ),
    left_annotation = rowAnnotation(
        `Source CT` = row_source_ct,
        col = list(`Source CT` = ct_colors),
        show_legend = FALSE,
        annotation_name_gp = gpar(fontsize = 8),
        simple_anno_size   = unit(0.4, "cm")
    ),
    right_annotation = rowAnnotation(
        label = anno_text(row_labels, gp = gpar(fontsize = 5.5), just = "left")
    ),
    width          = unit(22, "cm"),
    height         = unit(20, "cm"),
    use_raster     = TRUE,
    raster_quality = 2,
    heatmap_legend_param = list(title = "Z-score", at = c(-2, 0, 2), direction = "vertical")
)

pdf(file.path(plot_dir, "summary_heatmap.pdf"), width = 18, height = 14)
draw(hm2, merge_legend = TRUE, heatmap_legend_side = "right",
     annotation_legend_side = "right", padding = unit(c(4, 4, 4, 4), "mm"))
dev.off()

#-------------------------------------------------------------------------------
#   Correlation heatmap: peak-gene Pearson r for each trio x cell type
#
#   Each cell shows the correlation between a trio's peak accessibility and gene
#   expression across metacells of a given CT, using CT-Z-scored values.
#   Grey cells indicate zero-variance features in non-source CTs (NA).
#   Rows grouped by source CT and clustered within each group on all 3 features.
#   Called twice: once for all trios, once restricted to habenula neuron trios.
#   TF identity encoded as a row color annotation (top 20 TFs by count; rest gray).
#-------------------------------------------------------------------------------

# TF color palette: computed once from all trios so colors are consistent across
# both the full and Hb-subset versions of the heatmap
palette_20 <- c(
    "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4",
    "#42D4F4", "#F032E6", "#BFEF45", "#FABED4", "#469990",
    "#DCBEFF", "#9A6324", "#FFFAC8", "#800000", "#AAFFC3",
    "#808000", "#FFD8B1", "#000075", "#A9A9A9", "#000000"
)
# trio_sub: subset of trio_df; row_idx: corresponding row indices into feature matrices
make_cor_heatmap <- function(trio_sub, row_idx, filename) {
    source_ct_sub <- trio_sub$cell_type
    ct_sub_order  <- cell_type_order[cell_type_order %in% unique(source_ct_sub)]

    # Top 20 TFs and color palette computed from this specific subset
    top20_tfs <- trio_sub |> count(TF, sort = TRUE) |> slice_head(n = 20) |> pull(TF)
    tf_colors  <- c(setNames(palette_20, top20_tfs), "Other" = "gray80")

    # Row ordering: group by source CT, cluster within each group on all 3 features
    row_ord <- integer(0); row_grp <- character(0)
    for (ct in ct_sub_order) {
        ti  <- which(source_ct_sub == ct)
        idx <- ct_col_ranges[[ct]]
        sm  <- cbind(gene_mat_z2[row_idx[ti], idx, drop = FALSE],
                     peak_mat_z2[row_idx[ti], idx, drop = FALSE],
                     tf_mat_z2[row_idx[ti],   idx, drop = FALSE])
        o   <- if (length(ti) == 1) 1L else hclust(dist(sm), method = "complete")$order
        row_ord <- c(row_ord, ti[o])
        row_grp <- c(row_grp, rep(ct, length(ti)))
    }

    rsc    <- source_ct_sub[row_ord]
    tf_sub <- factor(
        ifelse(trio_sub$TF[row_ord] %in% top20_tfs, trio_sub$TF[row_ord], "Other"),
        levels = c(top20_tfs, "Other")
    )

    # Correlation matrix: rows in row_ord order, columns = all cell types
    # Restricted to metacells with nonzero raw TF expression (>= 3 required, else NA)
    cm <- matrix(NA_real_, nrow = nrow(trio_sub), ncol = length(ct_sub_order),
                 dimnames = list(NULL, ct_sub_order))
    for (ct in ct_sub_order) {
        idx <- ct_col_ranges[[ct]]
        for (i in seq_len(nrow(trio_sub))) {
            nz <- which(tf_mat[row_idx[row_ord[i]], idx] > 0)
            if (length(nz) < 3) next
            cm[i, ct] <- suppressWarnings(cor(
                gene_mat_z2[row_idx[row_ord[i]], idx[nz]],
                peak_mat_z2[row_idx[row_ord[i]], idx[nz]]
            ))
        }
    }

    hm <- Heatmap(
        cm,
        name   = "Peak-gene\ncorr.",
        col    = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        na_col = "grey85",
        cluster_rows       = FALSE,
        cluster_row_slices = FALSE,
        row_split          = factor(row_grp, levels = ct_sub_order),
        show_row_dend      = FALSE,
        show_row_names     = FALSE,
        row_title_gp       = gpar(fontsize = 7, fontface = "bold"),
        row_gap            = unit(1, "mm"),
        cluster_columns    = FALSE,
        show_column_names  = TRUE,
        column_names_gp    = gpar(fontsize = 7),
        column_names_rot   = 45,
        left_annotation  = rowAnnotation(
            `Source CT` = rsc,
            col = list(`Source CT` = ct_colors[ct_sub_order]),
            show_legend = FALSE,
            annotation_name_gp = gpar(fontsize = 8),
            simple_anno_size   = unit(0.4, "cm")
        ),
        right_annotation = rowAnnotation(
            TF  = tf_sub,
            col = list(TF = tf_colors),
            annotation_name_gp = gpar(fontsize = 8),
            simple_anno_size   = unit(0.4, "cm")
        ),
        width          = unit(8, "cm"),
        height         = unit(nrow(trio_sub) * 0.13, "cm"),
        use_raster     = TRUE,
        raster_quality = 2,
        heatmap_legend_param = list(title = "Peak-gene\ncorr.", at = c(-1, 0, 1), direction = "vertical")
    )

    pdf(file.path(plot_dir, filename), width = 10, height = 14)
    draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right",
         padding = unit(c(4, 4, 4, 4), "mm"))
    dev.off()
    invisible(NULL)
}

# All cell types
make_cor_heatmap(trio_df, seq_len(nrow(trio_df)), "peak_gene_correlation_heatmap.pdf")

# Habenula neurons only (LHb and MHb cell types)
hb_idx <- which(grepl("Hb", trio_df$cell_type))
make_cor_heatmap(trio_df[hb_idx, ], hb_idx, "peak_gene_correlation_heatmap_Hb.pdf")

session_info()
