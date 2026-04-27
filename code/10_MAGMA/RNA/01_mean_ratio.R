library(SingleCellExperiment)
library(tidyverse)
library(DeconvoBuddies)
library(here)
library(sessioninfo)
library(qs2)
library(rtracklayer)

sce_path = here(
    'processed-data', '05_03_annotation_adjustments', '06_refined_annotations',
    'refined_annotation_multiomeHab_SCE.qs2'
)
reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'
out_dir = here("processed-data", "10_MAGMA", "RNA", "gene_sets")
mean_ratio_threshold = 1.05
max_num_genes = 200

dir.create(out_dir, showWarnings = FALSE)

export_set = function(sce, cell_type_col, file_tag) {    
    marker_stats = get_mean_ratio(
            sce = sce, assay_name = "logcounts",
            cellType_col = cell_type_col, gene_ensembl = "gene_id",
            gene_name = "gene_name"
        ) |>
        filter(MeanRatio > mean_ratio_threshold) |>
        dplyr::rename(
            set_id = cellType.target, gene_id = gene_ensembl,
            mean_ratio = MeanRatio
        ) |>
        group_by(set_id) |>
        arrange(desc(mean_ratio)) |>
        slice_head(n = max_num_genes) |>
        ungroup() |>
        select(set_id, gene_id, gene_name, mean_ratio) |>
        arrange(set_id, desc(mean_ratio))

    message(sprintf("Marker counts for %s resolution:", file_tag))
    print(table(marker_stats$set_id))

    write_tsv(marker_stats, file.path(out_dir, sprintf("%s.tsv", file_tag)))
}

sce = qs_read(sce_path)

# add gene_id and gene_name to rowData
rowData(sce)$gene_name = rownames(sce)

gtf = import(reference_gtf)
gtf = gtf[gtf$type == "gene"]

rowData(sce)$gene_id <- gtf$gene_id[
    match(rowData(sce)$gene_name, gtf$gene_name)
]

# table(colData(sce)$refined_mid_cluster)

keep <- !is.na(rowData(sce)$gene_id) & rowData(sce)$gene_id != ""
message(sprintf("Dropping %d genes not in the GTF", sum(!keep)))
sce <- sce[keep, ]

sce$cell_type_mid = str_replace(
    sce$refined_mid_cluster, "(.*)([ML]Hb)(.*)", "\\2"
)
sce$cell_type_broad = str_replace(sce$cell_type_mid, "^[ML]Hb$", "Hb")

# export gene sets for both broad and fine resolutions
export_set(sce, "refined_mid_cluster", "fine")
export_set(sce, "cell_type_mid", "mid")
export_set(sce, "cell_type_broad", "broad")

session_info()
