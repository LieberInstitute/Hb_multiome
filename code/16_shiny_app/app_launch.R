library(Seurat)
library(qs2)
library(readr)

# For testing at JHPCE
# setwd(here::here("code", "16_shiny_app"))

source("app_core.R")

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
    Oligo = "#384a08",
    OPC = "#829454"
)

allowed_color_vars <- c("mid_cluster", "fine_cluster")

atlas_seur <- qs_read("atlas_seur_minimal.qs2")
metacell_seur <- qs_read("merged_metacell_seur.qs2")
trio_df <- read_csv("trios.csv.gz", show_col_types = FALSE)

atlas_seur@meta.data$mid_cluster <- factor(
    atlas_seur@meta.data$mid_cluster, levels = names(cell_type_colors)
)

run_app(
    atlas_seur = atlas_seur,
    color_vars = allowed_color_vars,
    metacell_seur = metacell_seur,
    trio_df = trio_df,
    cell_type_var = "mid_cluster",
    cell_type_colors = cell_type_colors,
    default_reduction = "wnn_umap"
)
