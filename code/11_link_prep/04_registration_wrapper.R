#   We'll at least need to register the Visium HD data to the latest multiome
#   (RNA) data, so compute the registration modeling statistics for multiome
#   here

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)
library(qs2)
library(rtracklayer)

sce_path = here(
    'processed-data', '05_03_annotation_adjustments', '06_refined_annotations',
    'refined_annotation_multiomeHab_SCE.qs2'
)
pseudo_path = here(
    'processed-data', '11_link_prep', '04_registration_wrapper',
    'pb_sce.rds'
)
model_path = here(
    'processed-data', '11_link_prep', '04_registration_wrapper',
    'modeling_results', 'model_results.rds'
)
reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'

dir.create(dirname(pseudo_path), showWarnings = FALSE)

sce = qs_read(sce_path)

#   Add in rowData with symbols and Ensembl IDs
gtf = import(reference_gtf) |>
    as.data.frame() |>
    filter(type == "gene") |>
    as_tibble() |>
    select(gene_id, gene_name)

rowData(sce)$gene_name = rownames(sce)
rowData(sce)$gene_id = gtf$gene_id[match(rowData(sce)$gene_name, gtf$gene_name)]
message(
    sprintf(
        "Removed %d genes without Ensembl IDs", sum(is.na(rowData(sce)$gene_id))
    )
)
sce = sce[!is.na(rowData(sce)$gene_id), ]
rownames(sce) = rowData(sce)$gene_id

#   Pseudobulk
model_results = registration_wrapper(
    sce,
    var_registration = 'refined_mid_cluster',
    var_sample_id = 'orig.ident',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path
)

saveRDS(model_results, model_path)

session_info()
