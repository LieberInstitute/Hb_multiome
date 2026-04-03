library(sessioninfo)
library(Seurat)
library(Signac)
library(tidyverse)
library(here)
library(qs2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(Matrix)

seur_in_path = here(
    'processed-data', '11_link_prep', '02_rebuild_atac_assay',
    'cell_level_seur.qs2'
)
pseudobulk_vars = c("refined_mid_cluster", "orig.ident")
prop_genes = 0.02

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
