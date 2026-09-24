library(here)
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)
library("SingleCellExperiment")
library("Seurat")
library(rtracklayer)
library(qs2)

#Seurat object for multiome

# RNAassey = here("processed-data", "05_Clustering_ARCr", "17_wnn_clustering_final_ct", "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds")
# rds<-readRDS(RNAassey)

# # change to SingleCellExperiment
# sce <- as.SingleCellExperiment(rds, assay = "RNA")

# *********************************************
sce <- qs_read('/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/05_03_annotation_adjustments/06_refined_annotations/refined_annotation_multiomeHab_SCE.qs2')

# registration wrapper

pseudo_path_broad = here(
    'processed-data', '10_MAGMA', 'RNA', 'registration_banksy',
    'pseudobulk_spe', 'broad.rds'
)
pseudo_path_fine = here(
    'processed-data', '10_MAGMA', 'RNA', 'registration_banksy',
    'pseudobulk_spe', 'fine.rds'
)

model_path_broad = here(
    'processed-data', '10_MAGMA', 'RNA', 'registration_banksy',
    'modeling_results', 'broad.rds'
)
model_path_fine = here(
    'processed-data', '10_MAGMA', 'RNA', 'registration_banksy',
    'modeling_results', 'fine.rds'
)

dir.create(dirname(pseudo_path_broad), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(model_path_broad), recursive = TRUE, showWarnings = FALSE)

rowData(sce)$gene_name = rownames(rowData(sce))

reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'
gtf = import(reference_gtf)
gtf = gtf[gtf$type == "gene"]

rowData(sce)$gene_id <- gtf$gene_id[match(rowData(sce)$gene_name, gtf$gene_name)]

table(colData(sce)$refined_mid_cluster)
table(colData(sce)$refined_cluster_ann)

keep <- !is.na(rowData(sce)$gene_id) & rowData(sce)$gene_id != ""
sce <- sce[keep, ]

model_results_broad = registration_wrapper(
    sce,
    var_registration = 'refined_mid_cluster',
    var_sample_id = 'orig.ident',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path_broad
)

model_results_fine = registration_wrapper(
    sce,
    var_registration = 'refined_cluster_ann',
    var_sample_id = 'orig.ident',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    pseudobulk_rds_file = pseudo_path_fine
)

saveRDS(model_results_broad, model_path_broad)
saveRDS(model_results_fine, model_path_fine)


