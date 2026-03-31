library(sessioninfo)
library(Seurat)
library(Signac)
library(GenomicRanges) 
library(tidyverse)
library(here)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(duckplyr)

seur_path = here(
    "processed-data", "06_peak_calling", "12_pseudobulk_MACS2",
    "Mid_pseudobulk.spearman.5e5.rds"
)
link_path = here(
    'processed-data', '12_new_peaks', '01_link_peaks', 'all_data.parquet'
)
plot_dir = here("plots", "12_new_peaks", "05_summary_heatmap")
atac_assay = "ATAC_macs2_pseudo"
rna_assay = "RNA"
highlight_cell_type = c("MHb.2", "LHb.2.7")
cor_thres = 0.3
FDR_thres = 0.1

cell_type_colors = c(
    "Oligo" = "#4d5802",
    "OPC" = "#d3c871",
    "Microglia" = "#222222",
    "Astrocyte" = "#8d363c",
    "Endo" = "#ee6c14",
    "Excit.Neuron" = "#666666",
    "Inhib.Thal" = '#b5a2ff',
    "Excit.Thal" = "#9e4ad1",
    "LHb.1.3.4" = "#0085af",
    "LHb.2.7" = "#7DF9FF",
    "LHb.4" = '#76c2af',
    "MHb.1" = "#FF00FF",
    "MHb.1.2" = "#c76a6a",
    "MHb.2" = "#FAA0A0",
    "MHb.3" = "#fa246a"
)

dir.create(plot_dir, showWarnings = FALSE)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

################################################################################
#   Functions
################################################################################

generate_summary_heatmap = function(
        count_df, cell_type_colors, plot_path, pdf_width = 10, pdf_height = 7
    ) {
    stopifnot(setequal(count_df$cell_type, names(cell_type_colors)))
    stopifnot(setequal(count_df$target_cell_type, names(cell_type_colors)))
  
    count_df = count_df |>
        group_by(row_id) |>
        mutate(
            peak_value = (peak_value - mean(peak_value)) / sd(peak_value),
            gene_value = (gene_value - mean(gene_value)) / sd(gene_value)
        ) |>
        ungroup() |>
        mutate(
            target_cell_type = factor(
                target_cell_type, levels = names(cell_type_colors)
            )
        ) |>
        arrange(target_cell_type)

    gene_mat_list = list()
    peak_mat_list = list()
    for (this_cell_type in names(cell_type_colors)) {
        gene_mat_list[[this_cell_type]] <- count_df |>
            filter(cell_type == this_cell_type) |>
            select(row_id, gene_value) |>
            distinct() |>
            column_to_rownames("row_id") |>
            as.matrix()

        peak_mat_list[[this_cell_type]] <- count_df |>
            filter(cell_type == this_cell_type) |>
            select(row_id, peak_value) |>
            distinct() |>
            column_to_rownames("row_id") |>
            as.matrix()
    }

    # Combine matrices: gene then peak for each cell type
    combined_mat <- do.call(
        cbind,
        lapply(
            names(cell_type_colors),
            function(ct) cbind(gene_mat_list[[ct]], peak_mat_list[[ct]])
        )
    )

    anno_df <- count_df |>
        select(row_id, target_cell_type) |>
        distinct()

    # Create color functions
    col_fun_gene <- colorRamp2(
        range(count_df$gene_value),
        c("#440154FF", "#FDE725FF")
    )

    col_fun_peak <- colorRamp2(
        range(count_df$peak_value),
        c("#000004FF", "#FCFDBFFF")
    )

    # Create column split and data type vectors
    cell_types <- names(cell_type_colors)
    col_split_labels <- unlist(lapply(cell_types, function(ct) {
        c(
            rep(paste0(ct, " - Gene"), ncol(gene_mat_list[[ct]])),
            rep(paste0(ct, " - Peak"), ncol(peak_mat_list[[ct]]))
        )
    }))
    col_split <- factor(
        col_split_labels,
        levels = unlist(lapply(
            cell_types, function(ct) paste0(ct, c(" - Gene", " - Peak"))
        ))
    )

    data_type_vec <- unlist(lapply(cell_types, function(ct) {
        c(
            rep("Gene", ncol(gene_mat_list[[ct]])),
            rep("Peak", ncol(peak_mat_list[[ct]]))
        )
    }))

    # Create column annotation
    cell_type_vec <- factor(
        unlist(lapply(cell_types, function(ct) {
            rep(ct, ncol(gene_mat_list[[ct]]) + ncol(peak_mat_list[[ct]]))
        })),
        levels = cell_types
    )

    data_type_colors <- c("Gene" = "#E69F00", "Peak" = "#56B4E9")

    col_ha <- HeatmapAnnotation(
        `Cell Type` = cell_type_vec,
        `Data Type` = data_type_vec,
        col = list(
            `Cell Type` = cell_type_colors,
            `Data Type` = data_type_colors
        ),
        show_legend = TRUE,
        show_annotation_name = FALSE
    )

    # Get row orderings by clustering on gene matrices per cell type
    full_row_order <- c()
    for (ct in cell_types) {
        rows_ct <- which(anno_df$target_cell_type == ct)
        if (length(rows_ct) > 0) {
            dummy_ht <- Heatmap(
                gene_mat_list[[ct]][rows_ct, , drop = FALSE],
                cluster_columns = FALSE
            )
            full_row_order <- c(
                full_row_order, rows_ct[row_order(dummy_ht)]
            )
        }
    }

    # Map colors to the combined matrix
    color_mat <- matrix(NA, nrow = nrow(combined_mat), ncol = ncol(combined_mat))
    for (i in seq_len(ncol(combined_mat))) {
        if (data_type_vec[i] == "Gene") {
            color_mat[, i] <- col_fun_gene(combined_mat[, i])
        } else {
            color_mat[, i] <- col_fun_peak(combined_mat[, i])
        }
    }

    # Create row annotation
    row_ha <- rowAnnotation(
        `Target Cell Type` = anno_df$target_cell_type,
        col = list(`Target Cell Type` = cell_type_colors),
        show_annotation_name = FALSE
    )

    # Create the main heatmap
    ht <- Heatmap(
        combined_mat,
        name = "Z-score",
        column_split = col_split,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_order = full_row_order,
        row_split = anno_df$target_cell_type,
        cluster_row_slices = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        row_title = "Link Cell Type",
        column_title = "Measured Cell Type",
        left_annotation = row_ha,
        top_annotation = col_ha,
        border = TRUE,
        show_heatmap_legend = FALSE,
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

    pdf(plot_path, width = pdf_width, height = pdf_height)
    draw(ht, annotation_legend_list = list(gene_legend, peak_legend))
    dev.off()
}

################################################################################
#   Main
################################################################################

seur = readRDS(seur_path)
seur@meta.data$donor = str_extract(colnames(seur), 'S[0-9]{2}.*')

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

# Calculate mean values across all cells for each cell type
count_df_list = list()
for (cell_type in unique(link_df$target_cell_type)) {
    subset_vec = (seur@meta.data$orig.ident == cell_type)
    
    count_df_list[[length(count_df_list) + 1]] <- tibble(
        cell_type = cell_type,
        target_cell_type = link_df$target_cell_type,
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

count_df = bind_rows(count_df_list) |>
    mutate(row_id = paste(peak, gene, sep = "|")) |>
    #   Temporary fix since the old object has extra cell types
    filter(
        target_cell_type %in% names(cell_type_colors),
        cell_type %in% names(cell_type_colors)
    )

cell_type_colors = cell_type_colors[
    names(cell_type_colors) %in% count_df$cell_type
]

#   All cell types
generate_summary_heatmap(
    count_df,
    cell_type_colors,
    plot_path = file.path(plot_dir, "all_cell_types.pdf"),
    pdf_width = 20,
    pdf_height = 18
)

#   Highlighted cell types
generate_summary_heatmap(
    count_df |>
        filter(
            target_cell_type %in% highlight_cell_type,
            cell_type %in% highlight_cell_type
        ),
    cell_type_colors[highlight_cell_type],
    plot_path = file.path(plot_dir, "highlighted_cell_types.pdf"),
    pdf_width = 10, pdf_height = 7
)

session_info()
  