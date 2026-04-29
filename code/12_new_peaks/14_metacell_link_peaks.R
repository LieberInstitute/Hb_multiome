#   As an array job over cell types, compute linked peaks using all peaks, even
#   those only called for other cell types

library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(Matrix)
library(sparseMatrixStats)
library(qs2)
library(duckplyr)

cell_types = c(
    "Astrocyte", "Endo", "Ependymal", "Excit.Thal", "Inhib_LHb_4.1",
    "Inhib_LHb_4.2", "Inhib.Thal", "LHb.1.3.4", "LHb.2.7", "LHb.4", "MHb.1",
    "MHb.1.2", "MHb.2", "MHb.3", "Microglia", "Oligo", "OPC"
)
this_cell_type = cell_types[as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))]

seur_path = here(
    "processed-data", "11_link_prep", "05_metacell_aggregate",
    "meta_seur.qs2"
)
out_path = here(
    "processed-data", "12_new_peaks", "14_metacell_link_peaks",
    sprintf("%s.parquet", this_cell_type)
)
min_cells = 6
min_prop_cells = 0.05
atac_assay = "ATAC"
rna_assay = "RNA"

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(out_path), showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

filter_features = function(seur, assay_name, min_cells) {
    counts_mat = GetAssayData(seur, assay = assay_name, layer = "counts")

    keep_features = rownames(counts_mat)[
        (Matrix::rowSums(counts_mat > 0) >= min_cells) &
        (sparseMatrixStats::rowSds(counts_mat) > 0)
    ]

    message(
        sprintf(
            "Filtering %s assay to %d of %d (%.1f%%) features (present in at least %d metacells)",
            assay_name, length(keep_features), nrow(counts_mat),
            100 * length(keep_features) / nrow(counts_mat), min_cells
        )
    )

    return(subset(seur[[assay_name]], features = keep_features))
}

################################################################################
#   Main
################################################################################

message(Sys.time(), ' | Loading Seurat object')
seur = qs_read(seur_path)

#   Subset to this cell type
seur = subset(seur, subset = refined_mid_cluster == this_cell_type)

#   Require both genes and peaks to be present in at least 6 metacells, or 5%,
#   whichever is larger
message(Sys.time(), ' | Filtering features')
this_min_cells = max(min_cells, as.integer(ncol(seur) * min_prop_cells))
seur[[atac_assay]] = filter_features(seur, atac_assay, this_min_cells)
seur[[rna_assay]] = filter_features(seur, rna_assay, this_min_cells)

#   Compute links without any filtering of outputs
message(Sys.time(), ' | Running LinkPeaks')
seur = LinkPeaks(
    object = seur, peak.assay = atac_assay, expression.assay = rna_assay,
    min.cells = this_min_cells, pvalue_cutoff = 1, score_cutoff = 0,
    method = "spearman"
)

#   Export linked peaks
message(Sys.time(), ' | Exporting to parquet')
Links(seur[[atac_assay]]) |>
    as.data.frame() |>
    as_tibble() |>
    mutate(
        cell_type = this_cell_type,
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
