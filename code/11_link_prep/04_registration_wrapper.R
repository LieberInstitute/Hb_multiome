#   We'll at least need to register the Visium HD data to the latest multiome
#   (RNA) data, so compute the registration modeling statistics for multiome
#   here

library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)
library(qs2)
library(rtracklayer)

resolution = c('mid', 'fine')[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))]

sce_path = here(
    'processed-data', '05_03_annotation_adjustments', '06_refined_annotations',
    'refined_annotation_multiomeHab_SCE.qs2'
)
pseudo_path = here(
    'processed-data', '11_link_prep', '04_registration_wrapper',
    sprintf('pb_sce_%s.rds', resolution)
)
model_path = here(
    'processed-data', '11_link_prep', '04_registration_wrapper',
    sprintf('model_results_%s.rds', resolution)
)
reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'

if (resolution == 'mid') {
    cluster_var = 'refined_mid_cluster'
} else {
    cluster_var = 'refined_cluster_ann'
}

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
    var_registration = cluster_var,
    var_sample_id = 'orig.ident',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path
)

saveRDS(model_results, model_path)

session_info()
