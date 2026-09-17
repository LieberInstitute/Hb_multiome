library(here)
library(tidyverse)
library(sessioninfo)

marker_script = here(
    "code", "04_DiffExpr_Clustering_seurat", "remote_DGE_marker_gene_lists.R"
)
out_path = here(
    "processed-data", "18_misc_supp_tables", "01_lit_data_driven_markers",
    "lit_data_driven_markers.csv"
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
source(marker_script)

a = get_Top50r_markers_genes_Hb()
dd_marker_df_list = list()
for (x in names(a)) {
    dd_marker_df_list[[x]] = tibble(
        gene = a[[x]], cell_type = x, marker_category = 'Data-driven'
    )
}

dd_marker_df = bind_rows(dd_marker_df_list) |>
    filter(grepl('[ML]Hb|Thal', cell_type)) |>
    mutate(cell_type = sub('^DD_', '', cell_type))

a = get_erik_and_Hb_markers_genes()
lit_marker_df_list = list()
for (x in names(a)) {
    lit_marker_df_list[[x]] = tibble(
        gene = a[[x]], cell_type = x, marker_category = 'Literature-based'
    )
}

lit_marker_df = bind_rows(lit_marker_df_list) |>
    filter(grepl('[ML]?H[Bb]|Thal', cell_type)) |>
    mutate(
        cell_type = str_extract(cell_type, '[ML]?H[Bb]|Thalamus') |>
            str_replace('HB', 'Hb')
    )

rbind(dd_marker_df, lit_marker_df) |>
    write_csv(out_path)

session_info()
