library(SingleCellExperiment)
library(tidyverse)
library(Seurat)
library(Signac)
library(here)
library(sessioninfo)
library(qs2)
library(lobstr)
library(scuttle)

seur_in_path = here(
    'processed-data', '16_shiny_app', '05_objects_for_sharing',
    'seur_cell_snMultiome_habenula_atlas.rds'
)
out_dir = here('processed-data', '16_shiny_app', 'cell_RNA_iSEE')

cell_type_colors = c(
    MHb_A = "#5e0c01",
    MHb_B = "#943f02", 
    MHb_C = "#f67104",
    MHb_D = "#f4d5ab",
    LHb_A = "#ee9630",
    LHb_B = "#306171",
    LHb_C = "#082844",
    GABA_LHb_C.1 = "#9c66c0",
    GABA_LHb_C.2 = "#5e0c56",
    Excit.Thal = "#2e6296",
    Inhib.Thal = "#8DADCA",
    Astrocyte = "#972f2f",
    Endo = "#f65a45",
    Ependymal = "#dbb369",
    Microglia = "#141b02",
    Oligo =  "#384a08",
    OPC = "#829454"
) 

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seur = readRDS(seur_in_path)

#   Build directly from the RNA data layer to avoid copying counts or ATAC.
rna_data = LayerData(seur, assay = 'RNA', layer = 'data')
cell_metadata = seur[[]]
cell_metadata$sample_id = seur$orig.ident
cell_metadata$orig.ident = NULL

sce = SingleCellExperiment(
    assays = list(logcounts = rna_data),
    rowData = S4Vectors::DataFrame(seur[['RNA']][[]]),
    colData = S4Vectors::DataFrame(cell_metadata)
)

sce$mid_cluster_colors = cell_type_colors[sce$mid_cluster]
sce$mid_cluster = factor(sce$mid_cluster, levels = names(cell_type_colors))
stopifnot(!any(is.na(sce$mid_cluster_colors)))

#   Match Seurat's SCE conversion, retaining all cell embeddings.
for (reduction_name in Reductions(seur)) {
    embedding = Embeddings(seur[[reduction_name]])
    stopifnot(identical(rownames(embedding), colnames(sce)))
    reducedDim(sce, reduction_name) = embedding
}

#   Use symbols for rownames but fall back on ENSEMBL for duplicates
rownames(sce) = uniquifyFeatureNames(
    rowData(sce)$gene_id, rowData(sce)$gene_name
)

print(obj_size(sce))
qs_save(sce, file.path(out_dir, 'sce_cell_RNA_logcounts.qs2'))

session_info()
