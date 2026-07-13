library(Seurat)
library(qs2)
library(readr)

source("code/16_shiny_app/app_core.R")

atlas_path <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/07_iSEE_app/01_prep_sce/seur_multiome_habenula_atlas.rds"
metacell_path <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/16_shiny_app/01_prep_objects/merged_metacell_seur.qs2"
trio_path <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/16_shiny_app/01_prep_objects/trios.csv"

cell_type_colors = c(
    Astrocyte = "#972f2f",
    OPC = "#829454",
    Oligo =  "#384a08",
    Microglia = "#141b02",
    Endo = "#f65a45",
    Ependymal = "#dbb369",
    
    Excit.Thal = "#2e6296",
    Inhib.Thal = "#8DADCA",
    LHb.4 = "#082844",
    LHb_C = "#082844",

    GABA_LHb_C.1 = "#9c66c0",
    GABA_LHb_C.2 = "#5e0c56",
    Inhib_LHb_4.1 = "#9c66c0",
    Inhib_LHb_4.2 = "#5e0c56",

    LHb.2.7 = "#ee9630",
    LHb_A = "#ee9630",

    LHb.1.3.4 = "#306171",
    LHb_B = "#306171",

    MHb.1 = "#5e0c01",
    MHb_A = "#5e0c01",

    MHb.1.2 = "#f67104",
    MHb_C = "#f67104",

    MHb.2 = "#943f02",
    MHb_B = "#943f02",
    
    MHb.3 = "#f4d5ab", 
    MHb_D = "#f4d5ab"
)

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
