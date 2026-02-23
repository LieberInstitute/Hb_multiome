#   To what degree is accessibility noticeable in cell types where peaks were
#   not called by MACS2?

library(sessioninfo)
library(Seurat)
library(Signac)
library(GenomicRanges) 
library(tidyverse)
library(here)
library(Matrix)

seur_path = here(
    "processed-data", "06_peak_calling", "12_pseudobulk_MACS2",
    "Mid_pseudobulk.spearman.5e5.rds"
)
atac_assay = "ATAC_macs2_pseudo"

seur = readRDS(seur_path)

#   Manipulate info about which peaks were called in which cell types, preparing
#   for a sparse matrix product for efficient computation

#   Get metadata with sample and cell type info
meta = seur@meta.data |>
    rownames_to_column("sample_id") |>
    mutate(sample_idx = row_number()) |>
    as_tibble() |>
    select(sample_id, orig.ident, sample_idx)

total_donors = length(unique(str_extract(meta$sample_id, '_S[0-9]{2}')))
total_cell_types = length(unique(meta$orig.ident))

#   Essentially create a tibble of indices used to populate a sparse indicator
#   matrix
peak_df = tibble(
        peak_id = rownames(seur[[atac_assay]]),
        peak_called_in = granges(seur[[atac_assay]])$peak_called_in
    ) |>
    mutate(
        cell_types = strsplit(peak_called_in, ",")
    ) |>
    unnest(cell_types) |>
    mutate(cell_type = str_trim(cell_types)) |>
    #   Here we leverage the many-to-many row-matching behavior to expand each
    #   combination of peak and sample for a given cell type
    left_join(
        meta,
        by = c("cell_type" = "orig.ident"),
        relationship = "many-to-many"
    ) |>
    #   But peaks are unique, so a simple match will work
    mutate(peak_idx = match(peak_id, rownames(seur[[atac_assay]]))) |>
    select(peak_id, cell_type, peak_idx, sample_idx)

#   Grab some scalars to later divide rowSums by-- will explain in later
#   comments. Since pseudbulking can drop donor x cell-type combinations,
#   scalars refer to the product of total donors x cell_types, counting raw
#   counts as 0 for missing combinations
count_scalar_df = peak_df |>
    select(peak_id, cell_type) |>
    distinct() |>
    group_by(peak_id) |>
    summarize(num_cell_types = n()) |>
    ungroup() |>
    slice(match(rownames(seur[[atac_assay]]), peak_id))
stopifnot(!any(is.na(count_scalar_df)))

target_count_scalar = count_scalar_df$num_cell_types * total_donors
other_count_scalar = total_cell_types * total_donors - target_count_scalar

#   Create sparse indicator matrix: peaks × (sample + cell type)s
indicator = sparseMatrix(
    i = peak_df$peak_idx,
    j = peak_df$sample_idx,
    x = 1,
    dims = c(nrow(seur[[atac_assay]]), ncol(seur))
)

#   Get sum of raw counts for peaks in cell types where they were called vs not
#   called. Since the underlying matrix operation is simply adding up raw counts
#   across the relevant samples, we have to divide by scalars computed from
#   above. The resulting quantity 'target_sum' is interpreted as the average
#   counts for a sample of the given cell type
atac_counts = GetAssayData(seur, assay = atac_assay, layer = "counts")
count_df = tibble(
        peak_id = rownames(seur[[atac_assay]]),
        total_counts = rowSums(atac_counts)
    ) |>
    mutate(
        target_mean = rowSums(atac_counts * indicator) / target_count_scalar,
        nontarget_mean = rowSums(atac_counts * (1 - indicator)) / other_count_scalar
    ) |>
    #   NAs can occur if a peak is called in every cell type (rare)
    mutate(nontarget_mean = ifelse(is.na(nontarget_mean), 0, nontarget_mean)) |>
    select(peak_id, target_mean, nontarget_mean)

message(
    sprintf(
        "Target mean higher in %.2f%% of peaks",
        mean(count_df$target_mean > count_df$nontarget_mean) * 100
    )
)
