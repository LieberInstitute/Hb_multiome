library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(qs2)
library(SingleCellExperiment)
library(rtracklayer)
library(lobstr)
library(scuttle)

seur_in_path = here(
    'processed-data', '11_link_prep', '02_rebuild_atac_assay',
    'cell_level_seur.qs2'
)
peak_path = here(
    'processed-data', '11_link_prep', '01_call_peaks',
    'macs3_peaks.csv.gz'
)
cell_map_path = here('raw-data', 'cell_type_map.csv')
gtf_path = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'
out_dir = here('processed-data', '07_iSEE_app', '01_prep_sce')
metadata_drop_cols = c(
    'high.tss', 'blacklist_fraction', 'blacklist_ratio', 'pct_reads_in_peaks',
    'cluster_ann', 'merged_cluster', 'mid_cluster', 'nucleosome_group'
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
#   Clean up dimensional reductions
################################################################################

seur = qs_read(seur_in_path)

message("Original Seurat object size: ")
print(obj_size(seur))

#   We have a huge number of ambiguously named dimensional reductions, and some
#   are on outdated data. Keep only harmonized reductions. Among these, keep the
#   raw reductions used for RNA and ATAC (PCA and LSI, respectively) and
#   UMAPs based on those reductions, as well as the WNN UMAP 
keep_map <- c(
    integrated.harmony     = "rna_pca_harmony",
    integrated.lsi.harmony = "atac_lsi_harmony",
    umap.integrated        = "rna_umap_harmony",
    umap.lsi.integrated    = "atac_umap_harmony",
    wnn.umap               = "wnn_umap"
)

for (i in seq_along(keep_map)) {
    this_reduction = seur[[names(keep_map)[i]]]
    Key(this_reduction) = paste0(str_replace_all(keep_map[i], '_', ''), "_")
    seur[[keep_map[i]]] = this_reduction
}

for (this_reduction_name in Reductions(seur)[!Reductions(seur) %in% keep_map]) {
    seur[[this_reduction_name]] = NULL
}

################################################################################
#   Clean up metadata and add peak_called_in column to ATAC assay
################################################################################

cluster_map = read_csv(cell_map_path, show_col_types = FALSE)
rename_map <- stats::setNames(
    cluster_map$new_cell_type, cluster_map$old_cell_type
)

seur@meta.data <- seur@meta.data |>
    select(-all_of(metadata_drop_cols)) |>
    dplyr::rename(
        fine_cluster = refined_cluster_ann, mid_cluster = refined_mid_cluster
    ) |>
    #   Map old cell-type names to new ones
    mutate(
        across(
            c(fine_cluster, mid_cluster),
            ~ dplyr::coalesce(unname(rename_map[.x]), .x)
        )
    )

#   Add 'peak_called_in' metadata column to the ATAC assay
seur[['ATAC']][[]]$peak_called_in = tibble(
        peak_id = rownames(seur[['ATAC']])
    ) |>
    left_join(
        read_csv(peak_path, show_col_types = FALSE) |>
            mutate(peak_id = sprintf("%s-%d-%d", seqnames, start, end)) |>
            select(peak_id, peak_called_in),
        by = "peak_id"
    ) |>
    pull(peak_called_in)

#   There are also genes with no expression, which we don't need
keep_genes <- rownames(seur[['RNA']])[
    rowSums(GetAssayData(seur, assay = "RNA", layer = "counts") > 0) > 0
]
seur[["RNA"]] <- subset(seur[["RNA"]], features = keep_genes)

################################################################################
#   Convert to individual SingleCellExperiments for each assay (for iSEE)
################################################################################

sce = as.SingleCellExperiment(seur, assay = 'RNA')
reducedDimNames(sce) = tolower(reducedDimNames(sce))
sce$ident = NULL

sce_atac = as.SingleCellExperiment(seur, assay = 'ATAC')
reducedDimNames(sce_atac) = tolower(reducedDimNames(sce_atac))
sce_atac$ident = NULL

assays(sce) = list(logcounts = assays(sce)$logcounts)
assays(sce_atac) = list(logcounts = assays(sce_atac)$logcounts)
gc()

################################################################################
#   Clean up rowData() and issues with gene names
################################################################################

#   For the iSEE app, we want ENSEMBL IDs, gene symbols, and additional
#   metadata from the GTF. If we naively join with the GTF by gene symbol, 9
#   genes fail to join. This is due to an upstream issue (that really should've
#   been fixed earlier) where because non-unique gene symbols are used for the
#   rownames, make.unique() automatically added ".1" and similar suffixes.

gtf = import(gtf_path) |>
    as.data.frame() |>
    as_tibble() |>
    filter(type == "gene")

#   Build ordered ENSEMBL ID lookup per gene symbol from GTF.
#   rownames(sce) may be corrupted by make.unique() which appends .1, .2, etc.
#   to duplicate gene symbols, creating nonexistent names. Legitimate gene names
#   that happen to end in .[0-9]+ are distinguished by their direct presence in
#   the GTF. Within each symbol group, the nth rowname is mapped to the nth
#   ENSEMBL ID seen for that symbol in the GTF.
gtf_symbol_ensembl <- gtf |>
    select(gene_name, gene_id) |>
    distinct()

gtf_all_symbols <- unique(gtf_symbol_ensembl$gene_name)

#   Named list: gene symbol → character vector of ENSEMBL IDs in GTF order
gtf_ensembls_by_symbol <- gtf_symbol_ensembl |>
    group_by(gene_name) |>
    summarise(gene_ids = list(gene_id), .groups = "drop") |>
    deframe()

original_rownames    <- rownames(sce)
times_symbol_seen    <- integer(0)   # named counter: ENSEMBL IDs assigned per symbol so far
resolved_symbols     <- character(length(original_rownames))
resolved_ensembl_ids <- character(length(original_rownames))

for (row_idx in seq_along(original_rownames)) {
    this_rn     <- original_rownames[row_idx]
    #   If the rowname is directly in the GTF it is a legitimate name (e.g.
    #   "AL627309.1"); otherwise strip a trailing make.unique() suffix to recover
    #   the base symbol (e.g. "TP53.1" → "TP53").
    base_symbol <- if (this_rn %in% gtf_all_symbols) this_rn else sub("\\.[0-9]+$", "", this_rn)
    stopifnot(base_symbol %in% gtf_all_symbols)

    prior_assign_count  <- times_symbol_seen[base_symbol]
    if (is.na(prior_assign_count)) prior_assign_count <- 0L
    ensembls_for_symbol <- gtf_ensembls_by_symbol[[base_symbol]]
    stopifnot(prior_assign_count + 1L <= length(ensembls_for_symbol))

    resolved_ensembl_ids[row_idx]  <- ensembls_for_symbol[prior_assign_count + 1L]
    resolved_symbols[row_idx]      <- base_symbol
    times_symbol_seen[base_symbol] <- prior_assign_count + 1L
}

#   Join full GTF metadata by gene_id (unique), avoiding the non-unique gene_name
rowData(sce) <- tibble(
        gene_id = resolved_ensembl_ids, gene_name = resolved_symbols,
    ) |>
    left_join(
        gtf |> distinct(gene_id, .keep_all = TRUE) |> select(-gene_name),
        by = "gene_id"
    ) |>
    select(where(~!all(is.na(.x)))) |>
    DataFrame()

rownames(sce) = uniquifyFeatureNames(
    rowData(sce)$gene_id, rowData(sce)$gene_name
)

#   Probably also makes sense to use this as the rownames for the Seurat object
#   as well (and rowData)
stopifnot(identical(original_rownames, rownames(seur[['RNA']])))
rownames(seur[['RNA']]) = rownames(sce)
seur[['RNA']]@meta.data = rowData(sce) |> as.data.frame()

################################################################################
#   Save objects
################################################################################

message("Final Seurat object size: ")
print(obj_size(seur))

message("iSEE object sizes (RNA then ATAC):")
print(obj_size(sce))
print(obj_size(sce_atac))

#   ExperimentHub
saveRDS(seur, file = file.path(out_dir, "seur_multiome_habenula_atlas.rds"))

#   iSEE
qs_save(sce, file = file.path(out_dir, "sce_RNA_iSEE.qs2"))
qs_save(sce_atac, file = file.path(out_dir, "sce_ATAC_iSEE.qs2"))

session_info()
