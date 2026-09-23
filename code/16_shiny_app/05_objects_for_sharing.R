library(here)
library(tidyverse)
library(qs2)
library(Seurat)
library(Signac)
library(sessioninfo)
library(rtracklayer)
library(scuttle)

cell_in_path = here(
    'processed-data', '11_link_prep', '02_rebuild_atac_assay',
    'cell_level_seur.qs2'
)
metacell_in_path = here(
    'processed-data', '16_shiny_app', '01_prep_objects',
    'merged_metacell_seur.qs2'
)
peak_path = here(
    'processed-data', '11_link_prep', '01_call_peaks',
    'macs3_peaks.csv.gz'
)
gtf_path = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'
out_dir = here('processed-data', '16_shiny_app', '05_objects_for_sharing')
metadata_drop_cols = c(
    'high.tss', 'blacklist_fraction', 'blacklist_ratio', 'pct_reads_in_peaks',
    'cluster_ann', 'merged_cluster', 'mid_cluster', 'nucleosome_group'
)
cell_map_path = here('raw-data', 'cell_type_map.csv')

dir.create(out_dir, showWarnings = FALSE)

################################################################################
#   Clean up dimensional reductions
################################################################################

seur = qs_read(cell_in_path)

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
    pull(peak_called_in) |>
    #   'peak_called_in' is a comma-separated list of old cell-type names;
    #   map each element to its new name before rejoining
    map_chr(
        ~ str_split(.x, ",")[[1]] |>
            (\(x) dplyr::coalesce(unname(rename_map[x]), x))() |>
            paste(collapse = ",")
    )

#   There are also genes with no expression, which we don't need
keep_genes <- rownames(seur[['RNA']])[
    rowSums(GetAssayData(seur, assay = "RNA", layer = "counts") > 0) > 0
]
seur[["RNA"]] <- subset(seur[["RNA"]], features = keep_genes)

################################################################################
#   Clean up (Seurat's equivalent of) rowData() and issues with gene names
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

original_rownames    <- rownames(seur[['RNA']])
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
new_row_data <- tibble(
        gene_id = resolved_ensembl_ids, gene_name = resolved_symbols,
    ) |>
    left_join(
        gtf |> distinct(gene_id, .keep_all = TRUE) |> select(-gene_name),
        by = "gene_id"
    ) |>
    select(where(~!all(is.na(.x)))) |>
    as.data.frame()
rownames(new_row_data) <- original_rownames

new_rownames = uniquifyFeatureNames(new_row_data$gene_id, new_row_data$gene_name)

#   Probably also makes sense to use this as the rownames for the Seurat object
#   as well (and rowData)
stopifnot(identical(original_rownames, rownames(seur[['RNA']])))
rownames(new_row_data) <- new_rownames
rownames(seur[['RNA']]) = new_rownames
seur[['RNA']]@meta.data = new_row_data

saveRDS(seur, file.path(out_dir, "seur_cell_snMultiome_habenula_atlas.rds"))

################################################################################
#   Clean up metacell Seurat object analogously
################################################################################

metacell_seur = qs_read(metacell_in_path)

#   Global metadata has the same old-cell-type-name issue, including in the
#   rownames themselves (of the form "<old_cell_type>_metacell_<N>")
old_prefixes <- str_remove(rownames(metacell_seur@meta.data), "_metacell_\\d+$")
metacell_suffixes <- str_extract(rownames(metacell_seur@meta.data), "_metacell_\\d+$")
stopifnot(!anyNA(old_prefixes), !anyNA(metacell_suffixes))
stopifnot(all(old_prefixes %in% names(rename_map)))

#   Renaming cells (rather than assigning rownames(meta.data) directly) keeps
#   the meta.data rownames and assay colnames in sync
metacell_seur <- RenameCells(
    metacell_seur,
    new.names = paste0(unname(rename_map[old_prefixes]), metacell_suffixes)
)
#   mid_cluster is a factor; recoding it directly with recode-style logic
#   would silently coerce it to character, so remap the levels instead
levels(metacell_seur@meta.data$mid_cluster) <- dplyr::coalesce(
    unname(rename_map[levels(metacell_seur@meta.data$mid_cluster)]),
    levels(metacell_seur@meta.data$mid_cluster)
)

#   'orig.ident' has the same old-name issue, but is additionally coarser
#   than the cell-level fine-grained grouping (it collapses
#   Inhib_LHb_4.1/Inhib_LHb_4.2 into "Inhib"). The rowname-derived prefix is
#   the correct, fine-grained source of truth
metacell_seur@meta.data <- metacell_seur@meta.data |>
    mutate(orig.ident = factor(unname(rename_map[old_prefixes])))

#   rowData for the RNA assay is currently empty; populate it to match the
#   cell-level object, since the metacell RNA features are a subset of the
#   cell-level RNA features. Metacell RNA rownames still use the *original*
#   (pre-fix) gene symbols, so match against those, then relabel with the
#   final, uniquified gene names used in the cell-level object
match_idx_rna <- match(rownames(metacell_seur[['RNA']]), original_rownames)
stopifnot(!anyNA(match_idx_rna))

rownames(metacell_seur[['RNA']]) <- new_rownames[match_idx_rna]
metacell_seur[['RNA']][[]] <- new_row_data[match_idx_rna, ]

#   Add 'peak_called_in' metadata column to the ATAC assay, as with the
#   cell-level object
metacell_seur[['ATAC']][[]]$peak_called_in <- tibble(
        peak_id = rownames(metacell_seur[['ATAC']])
    ) |>
    left_join(
        read_csv(peak_path, show_col_types = FALSE) |>
            mutate(peak_id = sprintf("%s-%d-%d", seqnames, start, end)) |>
            select(peak_id, peak_called_in),
        by = "peak_id"
    ) |>
    pull(peak_called_in) |>
    map_chr(
        ~ str_split(.x, ",")[[1]] |>
            (\(x) dplyr::coalesce(unname(rename_map[x]), x))() |>
            paste(collapse = ",")
    )

saveRDS(
    metacell_seur,
    file.path(out_dir, "seur_metacell_snMultiome_habenula_atlas.rds")
)

session_info()
