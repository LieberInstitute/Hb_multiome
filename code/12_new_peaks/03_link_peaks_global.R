library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(Matrix)
library(sparseMatrixStats)
library(qs2)
library(duckplyr)

seur_path = here(
    'processed-data', '11_link_prep', '03_pseudobulk', 'pb_seur.qs2'
)
out_path = here(
    "processed-data", "12_new_peaks", "01_link_peaks",
    "global.parquet"
)
min_cells_prop = 0.1
atac_assay = "ATAC"
rna_assay = "RNA"

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

################################################################################
#   Functions
################################################################################

filter_features = function(seur, assay_name, min_cells_prop) {
    counts_mat = GetAssayData(seur, assay = assay_name, layer = "counts")
    keep_features = rownames(counts_mat)[
        (Matrix::rowSums(counts_mat > 0) >= min_cells_prop * ncol(counts_mat)) &
        (sparseMatrixStats::rowSds(counts_mat) > 0)
    ]

    message(
        sprintf(
            "Filtering %s assay to %d of %d (%.1f%%) features (present in at least %.1f%% of pb cells)",
            assay_name, length(keep_features), nrow(counts_mat),
            100 * length(keep_features) / nrow(counts_mat), min_cells_prop * 100
        )
    )

    return(subset(seur[[assay_name]], features = keep_features))
}

################################################################################
#   Main
################################################################################

message(Sys.time(), ' | Loading Seurat object')
seur = qs_read(seur_path)

#   Require both genes and peaks to be present in 10% of pseudobulked cells
message(Sys.time(), ' | Filtering features')
seur[[atac_assay]] = filter_features(seur, atac_assay, min_cells_prop)
seur[[rna_assay]] = filter_features(seur, rna_assay, min_cells_prop)

#   Compute links without any filtering of outputs
message(Sys.time(), ' | Running LinkPeaks')
seur = LinkPeaks(
    object = seur, peak.assay = atac_assay, expression.assay = rna_assay,
    min.cells = as.integer(ncol(seur) * min_cells_prop), pvalue_cutoff = 1,
    score_cutoff = 0, method = "spearman"
)

#   Export linked peaks
message(Sys.time(), ' | Exporting to parquet')
Links(seur[[atac_assay]]) |>
    as.data.frame() |>
    as_tibble() |>
    mutate(
        FDR = p.adjust(pvalue, method = "BH"),
        #  Avoids duckplyr issues
        across(where(is.factor), as.character)
    ) |>
    compute_parquet(out_path)

message("Memory usage:")
gc()

session_info()

## This script was made using slurmjobs version 1.3.0
## available from http://research.libd.org/slurmjobs/
