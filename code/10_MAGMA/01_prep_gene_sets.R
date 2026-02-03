#   Prepare gene sets for MAGMA. Each combination of cell type and peak category
#   will be a separate gene set (e.g. astrocyte enriched linked DAR)

library(tidyverse)
library(here)
library(sessioninfo)

pair_path = here(
    'processed-data', '06_peak_calling', '21_overlaping_FDRscores_TopHeatmap',
    'all_peaks_categorized.csv.gz'
)
out_path = here('processed-data', '10_MAGMA', 'gene_sets', '%s.tsv')

dir.create(dirname(out_path), showWarnings = FALSE)

pair_df = read_csv(pair_path, show_col_types = FALSE) |>
    #   Basically filter out unlinked DARs
    filter(!is.na(link_gene_id)) |>
    mutate(
        set_id = paste(
            cell_type,
            category |>
                str_replace_all(' +', '_'),
            sep = '_'
        )
    ) |>
    select(set_id, link_gene_id) |>
    #   At this point it's still posible for genes to be duplicated within a set
    #   (peaks can be different)
    distinct()

#   Gene sets at mid resolution
pair_df |>
    arrange(set_id) |>
    write_tsv(sprintf(out_path, 'mid'))

#   Gene sets at semi-broad resolution
pair_df |>
    mutate(set_id = str_replace(set_id, '^([ML])Hb\\.[^_]+', '\\1Hb')) |>
    distinct() |>
    arrange(set_id) |>
    write_tsv(sprintf(out_path, 'semi_broad'))

#   Gene sets at broad resolution
pair_df |>
    mutate(set_id = str_replace(set_id, '^[ML]Hb\\.[^_]+', 'Hb')) |>
    distinct() |>
    arrange(set_id) |>
    write_tsv(sprintf(out_path, 'broad'))

session_info()
