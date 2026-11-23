library(here)
library(tidyverse)
library(sessioninfo)

marker_path = here(
    'processed-data', '10_MAGMA', 'RNA', 'gene_sets', 'fine.tsv'
)
map_path = here('raw-data', 'cell_type_map.csv')

cell_type_colors = c(
    MHb_A = '#5e0c01',
    MHb_B = '#943f02',
    MHb_C = '#f67104',
    MHb_D = '#f4d5ab',
    LHb_A = '#ee9630',
    LHb_B = '#306171',
    LHb_C = '#082844',
    GABA_LHb_C.1 = '#9c66c0',
    GABA_LHb_C.2 = '#5e0c56',
    Excit.Thal = '#2e6296',
    Inhib.Thal = '#8DADCA',
    Astrocyte = '#972f2f',
    Endo = '#f65a45',
    Ependymal = '#dbb369',
    Microglia = '#141b02',
    Oligo = '#384a08',
    OPC = '#829454'
)

markers = read_tsv(marker_path, show_col_types = FALSE)
cell_map = read_csv(map_path, show_col_types = FALSE)

#   Require a one-to-one mapping between marker sets and the app's cell types.
stopifnot(
    setequal(unique(markers$set_id), cell_map$old_cell_type),
    setequal(cell_map$new_cell_type, names(cell_type_colors))
)

top_markers = markers |>
    left_join(cell_map, by = c('set_id' = 'old_cell_type')) |>
    mutate(mid_cluster = factor(new_cell_type, levels = names(cell_type_colors))) |>
    group_by(mid_cluster) |>
    arrange(desc(mean_ratio), gene_id, .by_group = TRUE) |>
    slice_head(n = 1) |>
    ungroup() |>
    arrange(mid_cluster)

message('Top markers by cell type:')
top_markers |>
    select(mid_cluster, gene_id, gene_name, mean_ratio) |>
    print(n = Inf)
top_markers |> pull(gene_name) |> dput()

session_info()
