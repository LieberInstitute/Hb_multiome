library(here)
library(tidyverse)
library(qs2)
library(Seurat)
library(Signac)
library(sessioninfo)
library(Matrix)
library(lobstr)

in_path = here(
    "processed-data", "07_iSEE_app", "01_prep_sce",
    "seur_multiome_habenula_atlas.rds"
)
out_dir = here("processed-data", "16_shiny_app", "02_prep_cell_level")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
#   Subset to a minimal Seurat with just reduced dims and metadata
################################################################################

seur = readRDS(in_path)

message("Original Seurat object size: ")
print(obj_size(seur))

empty_counts = Matrix(
    0,
    nrow = 1,
    ncol = ncol(seur),
    sparse = TRUE,
    dimnames = list("placeholder", colnames(seur))
)
placeholder_assay = CreateAssayObject(counts = empty_counts)
Key(placeholder_assay) = "placeholder_"

seur@assays = list(placeholder = placeholder_assay)
seur@active.assay = "placeholder"
seur@graphs = list()
seur@neighbors = list()
seur@commands = list()
seur@tools = list()
seur@misc = list()

for (nm in names(seur@reductions)) {
    seur@reductions[[nm]]@assay.used = "placeholder"
}

message("Final Seurat object size: ")
print(obj_size(seur))

qs_save(seur, file.path(out_dir, "atlas_seur_minimal.qs2"))

session_info()
