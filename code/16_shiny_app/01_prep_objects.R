library(tidyverse)
library(Seurat)
library(Signac)
library(qs2)
library(here)
library(sessioninfo)
library(Matrix)

cell_types = c(
    "Astrocyte", "Ependymal", "Excit.Thal", "Inhib_LHb_4.1",
    "Inhib_LHb_4.2", "Inhib.Thal", "LHb.1.3.4", "LHb.2.7", "LHb.4", "MHb.1",
    "MHb.1.2", "MHb.2", "Microglia", "Oligo", "OPC"
)
in_paths = here(
    "processed-data", "13_tripod_trios", "02_tripod_preprocess",
    sprintf("preprocessed_objects_%s.qs2", cell_types)
)
out_path = here(
    "processed-data", "16_shiny_app", "01_prep_objects",
    "merged_metacell_seur.qs2"
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

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

meta_df = tibble(
    cell = rownames(rna_mat),
    broad_cell_type = stringr::str_remove(cell, "_metacell_\\d+$")
) |>
    tibble::column_to_rownames("cell")

stopifnot(identical(rownames(rna_mat), rownames(peak_mat)))

merged_metacell_seur = CreateSeuratObject(
    counts = t(rna_mat),
    meta.data = meta_df,
    assay = "RNA"
)
merged_metacell_seur[["ATAC"]] = CreateAssayObject(counts = t(peak_mat))
DefaultAssay(merged_metacell_seur) = "RNA"

qs_save(merged_metacell_seur, out_path)

session_info()
