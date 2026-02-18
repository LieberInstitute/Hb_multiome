library(SingleCellExperiment)
library(tidyverse)
library(DeconvoBuddies)
library(here)
library(sessioninfo)

sce_path = here("processed-data", "10_MAGMA", "RNA", "sce.rds")
out_dir = here("processed-data", "10_MAGMA", "RNA", "gene_sets")
mean_ratio_threshold = 1.05
max_num_genes = 200

dir.create(out_dir, showWarnings = FALSE)

export_set = function(spe, cell_type_col, file_tag) {    
    marker_stats = get_mean_ratio(
            sce = spe_pb, assay_name = "logcounts",
            cellType_col = cell_type_col, gene_ensembl = "gene_id",
            gene_name = "gene_name"
        ) |>
        filter(MeanRatio > mean_ratio_threshold) |>
        dplyr::rename(set_id = cellType.target, gene_id = gene) |>
        group_by(set_id) |>
        arrange(desc(MeanRatio)) |>
        slice_head(n = max_num_genes) |>
        select(set_id, gene_id) |>
        arrange(set_id)

    message(sprintf("Marker counts for %s resolution:", file_tag))
    print(table(marker_stats$set_id))

    write_tsv(marker_stats, file.path(out_dir, sprintf("%s.tsv", file_tag)))
}

sce = readRDS(sce_path)

export_set(sce, "cell_type_fine", "fine")
export_set(sce, "cell_type_mid", "mid")
export_set(sce, "cell_type_broad", "broad")

session_info()
