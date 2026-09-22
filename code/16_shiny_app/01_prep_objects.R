library(tidyverse)
library(Seurat)
library(Signac)
library(qs2)
library(here)
library(sessioninfo)
library(Matrix)
library(duckplyr)

cell_types = c(
    "Astrocyte", "Ependymal", "Excit.Thal", "Inhib_LHb_4.1",
    "Inhib_LHb_4.2", "Inhib.Thal", "LHb.1.3.4", "LHb.2.7", "LHb.4", "MHb.1",
    "MHb.1.2", "MHb.2", "Microglia", "Oligo", "OPC"
)
cell_type_levels = c(
    'MHb_A', 'MHb_B', 'MHb_C', 'MHb_D', 'LHb_A', 'LHb_B', 'LHb_C',
    'GABA_LHb_C.1', 'GABA_LHb_C.2', 'Excit.Thal', 'Inhib.Thal', 'Astrocyte',
    'Endo', 'Ependymal', 'Microglia', 'Oligo', 'OPC'
)
dar_path = here("processed-data", "15_DARs", "03_gather", "DARs_fine.csv.gz")
cell_map_path = here('raw-data', 'cell_type_map.csv')
trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
in_paths = here(
    "processed-data", "13_tripod_trios", "02_tripod_preprocess",
    sprintf("preprocessed_objects_%s.qs2", cell_types)
)
out_dir = here("processed-data", "16_shiny_app", "01_prep_objects")

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
#   Prep metacell Seurat object
################################################################################

metacell_objs = purrr::map(in_paths, \(in_path) {
    qs_read(in_path)$metacell_seur
})
names(metacell_objs) = cell_types

merge_metacell_assay = function(metacell_objs, assay_name) {
    assay_mats = purrr::imap(metacell_objs, \(obj, cell_type) {
        mat = obj[[assay_name]]

        if (is.null(rownames(mat)) || is.null(colnames(mat))) {
            stop(sprintf("%s assay in %s is missing row or column names", assay_name, cell_type))
        }

        rownames(mat) = paste(cell_type, rownames(mat), sep = "_")
        Matrix(mat, sparse = TRUE)
    })

    all_features = assay_mats |>
        purrr::map(colnames) |>
        unlist(use.names = FALSE) |>
        unique()

    padded_mats = purrr::map(assay_mats, \(mat) {
        missing_features = setdiff(all_features, colnames(mat))

        if (length(missing_features) > 0) {
            zero_block = Matrix(
                0,
                nrow = nrow(mat),
                ncol = length(missing_features),
                dimnames = list(rownames(mat), missing_features),
                sparse = TRUE
            )
            mat = cbind(mat, zero_block)
        }

        mat[, all_features, drop = FALSE]
    })

    do.call(rbind, padded_mats)
}

rna_mat = merge_metacell_assay(metacell_objs, "rna")
peak_mat = merge_metacell_assay(metacell_objs, "peak")

cluster_map = read_csv(cell_map_path, show_col_types = FALSE)
rename_map <- stats::setNames(
    cluster_map$new_cell_type, cluster_map$old_cell_type
)

meta_df = tibble(
        cell = rownames(rna_mat),
        mid_cluster = stringr::str_remove(cell, "_metacell_\\d+$")
    ) |>
    tibble::column_to_rownames("cell") |>
    mutate(
        mid_cluster = dplyr::coalesce(
            unname(rename_map[mid_cluster]), mid_cluster
        ) |>
        factor(levels = cell_type_levels)
    )

stopifnot(identical(rownames(rna_mat), rownames(peak_mat)))

#   Metacell matrices take raw counts, library-size normalize them, multiply
#   by an arbitrary constant, but don't log-transform. Log-transform here.
#   https://github.com/yuchaojiang/TRIPOD/blob/c8421358f6bfddeb6b219cbc836a40aae9d20016/package/R/prep.R#L167-L171
rna_data  <- log1p(t(rna_mat))
atac_data <- log1p(t(peak_mat))

rna_assay  <- CreateAssay5Object(data = rna_data)
atac_assay <- CreateAssay5Object(data = atac_data)

seur <- CreateSeuratObject(
    counts = Matrix(
        0, nrow = nrow(rna_data), ncol = ncol(rna_data), sparse = TRUE,
        dimnames = dimnames(rna_data)
    ),
    meta.data = meta_df,
    assay = "RNA"
)

seur[["RNA"]] <- rna_assay
seur[["ATAC"]] <- atac_assay
DefaultAssay(seur) <- "RNA"

#   Add info about DARs to ATAC metadata
dar_df = read_csv(dar_path, show_col_types = FALSE) |>
    dplyr::rename(
        DAR_cell_type = cell_type, DAR_logFC = avg_log2FC, DAR_p_adj = p_val_adj
    ) |>
    select(peak, DAR_cell_type, DAR_logFC, DAR_p_adj)

seur[['ATAC']][[]] = tibble(
        peak = rownames(seur[['ATAC']]),
        peak_order = seq_along(rownames(seur[['ATAC']]))
    ) |>
    left_join(dar_df, by = 'peak') |>
    arrange(peak_order) |>
    select(-peak_order) |>
    column_to_rownames('peak')

qs_save(seur, file.path(out_dir, "merged_metacell_seur.qs2"))

################################################################################
#   Prep trios
################################################################################

trio_df = read_parquet_duckdb(trio_path, prudence = "stingy") |>
    filter(!is_intersect | (is_intersect & stringency_level == 1)) |>
    dplyr::rename(
        trio_coef = coef, trio_p_adj = adj, mid_cluster = cell_type,
        trio_level = stringency_level
    ) |>
    collect() |>
    mutate(
        trio_level = factor(
            ifelse(is_intersect, 'Intersect', trio_level),
            levels = c('Intersect', '1', '2')
        ),
        mid_cluster = factor(
            dplyr::coalesce(unname(rename_map[mid_cluster]), mid_cluster),
            levels = cell_type_levels
        )
    ) |>
    select(trio_level, peak, gene, TF, mid_cluster, trio_coef, trio_p_adj) |>
    arrange(trio_level, mid_cluster, trio_p_adj, trio_coef) |>
    write_csv(file.path(out_dir, "trios.csv.gz"))

#   Average (log-normalized) RNA expression per gene, within each cell type,
#   for annotating the gene and TF in each trio
rna_data_mat = seur[['RNA']]$data
avg_rna_expr = sapply(levels(seur$mid_cluster), function(ct) {
    cells = rownames(seur@meta.data)[seur$mid_cluster == ct]
    if (length(cells) == 0) return(rep(NA_real_, nrow(rna_data_mat)))
    Matrix::rowMeans(rna_data_mat[, cells, drop = FALSE])
})
rownames(avg_rna_expr) = rownames(rna_data_mat)

get_avg_expr = function(gene, mid_cluster) {
    row_idx = match(gene, rownames(avg_rna_expr))
    col_idx = match(as.character(mid_cluster), colnames(avg_rna_expr))
    ifelse(
        is.na(row_idx) | is.na(col_idx), NA_real_,
        avg_rna_expr[cbind(row_idx, col_idx)]
    )
}

#   Alternative export: wide form preserving both level-1 and level-2 stats.
#   Intersect trios have all four columns populated; non-intersect trios have
#   NAs in stats for one of the trio levels. This version is for a supplemental
#   table
read_parquet_duckdb(trio_path, prudence = "stingy") |>
    dplyr::rename(
        trio_coef = coef, trio_p_adj = adj, mid_cluster = cell_type,
        trio_level = stringency_level
    ) |>
    collect() |>
    mutate(
        mid_cluster = factor(
            dplyr::coalesce(unname(rename_map[mid_cluster]), mid_cluster),
            levels = cell_type_levels
        )
    ) |>
    pivot_wider(
        id_cols = c(peak, gene, TF, mid_cluster, is_intersect, is_unique_intersect, is_top_TF),
        names_from = trio_level,
        values_from = c(trio_coef, trio_p_adj),
        names_glue = "{.value}_lvl{trio_level}"
    ) |>
    mutate(
        gene_avg_expr = get_avg_expr(gene, mid_cluster),
        TF_avg_expr = get_avg_expr(TF, mid_cluster)
    ) |>
    arrange(mid_cluster, trio_p_adj_lvl1, trio_coef_lvl1) |>
    select(-c(is_intersect, is_unique_intersect, is_top_TF)) |>
    write_csv(file.path(out_dir, "trios_supp_table.csv.gz"))


session_info()
