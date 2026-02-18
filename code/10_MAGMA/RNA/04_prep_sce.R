#   Convert the huge multiome Seurat object into a smaller SingleCellExperiment
#   with Ensembl gene IDs

library(here)
library(sessioninfo)
library(tidyverse)
library(SingleCellExperiment)
library(Seurat)
library(Signac)
library(rtracklayer)

seurat_path = here(
    "processed-data", "05_Clustering_ARCr", "17_wnn_clustering_final_ct",
    "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
)
out_path = here("processed-data", "10_MAGMA", "RNA", "sce.rds")
reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'

#   Load and convert to SingleCellExperiment
sce = as.SingleCellExperiment(readRDS(seurat_path), assay = "RNA")

#   Trim object
assays(sce)$scaledata = NULL
reducedDims(sce) = list()

#   Remove genes with no expression
sce = sce[rowSums(assays(sce)$counts) > 0, ]

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

#   Define cell types at different resolutions
sce$cell_type_fine = sub("^C\\.\\d+\\.", "", sce$ident)
sce$cell_type_mid = str_replace(
    sce$cell_type_fine, '^([ML])Hb\\.[^_]+', '\\1Hb'
)
sce$cell_type_broad = str_replace(sce$cell_type_fine, '^[ML]Hb\\.[^_]+', 'Hb')

saveRDS(sce, out_path)

session_info()
