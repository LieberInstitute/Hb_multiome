library(Seurat)
library(qs2)
library(readr)

source("code/16_shiny_app/app_core.R")

atlas_path <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/07_iSEE_app/01_prep_sce/seur_multiome_habenula_atlas.rds"
metacell_path <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/16_shiny_app/01_prep_objects/merged_metacell_seur.qs2"
trio_path <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/16_shiny_app/01_prep_objects/trios.csv"

allowed_color_vars <- c("mid_cluster", "fine_cluster")

atlas_seur <- readRDS(atlas_path)
metacell_seur <- qs_read(metacell_path)
trio_df <- read_csv(trio_path, show_col_types = FALSE)

run_app(
    atlas_seur = atlas_seur,
    color_vars = allowed_color_vars,
    metacell_seur = metacell_seur,
    trio_df = trio_df,
    default_reduction = "wnn_umap",
    default_color_by = "mid_cluster"
)
