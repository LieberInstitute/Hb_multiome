library(Seurat)

source("code/16_shiny_app/app_core.R")

seurat_path <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/07_iSEE_app/01_prep_sce/seur_multiome_habenula_atlas.rds"
allowed_color_vars <- c(
  "orig.ident",
  "scDblFinder.class",
  "mid_cluster",
  "fine_cluster"
)

seur <- readRDS(seurat_path)

run_app(
  seur = seur,
  color_vars = allowed_color_vars,
  default_reduction = "wnn_umap"
)
