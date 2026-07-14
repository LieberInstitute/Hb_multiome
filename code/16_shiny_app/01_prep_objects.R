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
        )
    )

stopifnot(identical(rownames(rna_mat), rownames(peak_mat)))

seur = CreateSeuratObject(
    counts = t(rna_mat),
    meta.data = meta_df,
    assay = "RNA"
)
seur[["ATAC"]] = CreateAssayObject(counts = t(peak_mat))
DefaultAssay(seur) = "RNA"

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
    filter(stringency_level == 1, is_intersect) |>
    dplyr::rename(trio_cor = coef, trio_p_adj = adj) |>
    select(peak, gene, TF, cell_type, trio_cor, trio_p_adj) |>
    collect() |>
    write_csv(file.path(out_dir, "trios.csv"))

session_info()
