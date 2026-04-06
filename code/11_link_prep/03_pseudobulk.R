library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(qs2)
library(Matrix)
library(edgeR)

seur_in_path = here(
    'processed-data', '11_link_prep', '02_rebuild_atac_assay',
    'cell_level_seur.qs2'
)
seur_out_path = here(
    'processed-data', '11_link_prep', '03_pseudobulk', 'pb_seur.qs2'
)
pseudobulk_vars = c("refined_mid_cluster", "orig.ident")
prop_genes = 0.02

dir.create(dirname(seur_out_path), showWarnings = FALSE)

seur = qs_read(seur_in_path)

#   As Cynthia originally did, filter genes to those in a minimum proportion of cells
keep_genes = rownames(seur[["RNA"]])[
    Matrix::rowSums(seur[["RNA"]]@counts > 0) > prop_genes * ncol(seur[["RNA"]])
]
seur[["RNA"]] = subset(seur[["RNA"]], features = keep_genes)

seur_pb = AggregateExpression(
    seur, assays = c("RNA", "ATAC"),
    group.by = pseudobulk_vars,
    return.seurat = FALSE,
    verbose = TRUE
)

#   Use bulk-style normalization (logCPM) for ATAC, since the TFIDF
LayerData(seur_pb, assay = "ATAC", layer = "data") = GetAssayData(
        seur_pb, assay = "ATAC", layer = "counts"
    ) |>
    calcNormFactors() |>
    edgeR::cpm(log = TRUE, prior.count = 1)

seur_pb[["RNA"]] = NormalizeData(seur_pb[["RNA"]])

qs_save(seur_pb, seur_out_path)

session_info()
