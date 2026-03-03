#   Write an all-purpose Seurat object to disk with some data trimmed, and using
#   serialization with qs2. The main priority is lowering load time for the
#   Seurat object, which takes ~30min (!) before this script

library(Seurat)
library(Signac)
library(here)
library(qs2)
library(sessioninfo)

in_path = here(
    "processed-data", "05_Clustering_ARCr", "22_add_mid_level_clustering",
    "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
)
out_path = here(
    "processed-data", "05_Clustering_ARCr", "24_smaller_seurat",
    "general_purpose_seur.qs2"
)
keep_reductions = c("pca", "umap.integrated", "umap.lsi.integrated", "wnn.umap")

message(Sys.time(), " | Reading original object")
seur = readRDS(in_path)

#   Drop most reductions
message(Sys.time(), " | Dropping some data")
remove_reductions = setdiff(names(seur@reductions), keep_reductions)
for (red in remove_reductions) {
    seur[[red]] = NULL
}

#   Drop 'scale.data' from the RNA assay
seur[["RNA"]] = subset(
    seur[["RNA"]], cells = colnames(seur), layer = c("data", "counts")
)

message(Sys.time(), " | Saving smaller Seurat object")
qs_save(seur, out_path)

session_info()
