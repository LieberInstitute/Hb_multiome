#   Prepare gene sets for MAGMA. Each combination of cell type and peak category
#   will be a separate gene set (e.g. astrocyte enriched linked DAR)

library(tidyverse)
library(here)
library(sessioninfo)

pair_path = here(
    'processed-data', '06_peak_calling', '21_overlaping_FDRscores_TopHeatmap',
    'overlaps_linkPeak_DARs_classified_thr_CC0.3_thr_DAR0.1.csv'
)
out_path = here(
    'processed-data', '10_MAGMA', 'first_test', 'input_gene_sets.tsv'
)

read_csv(pair_path, show_col_types = FALSE) |>
    filter(
        type_classification == 'classification-3', category != "Non-significant"
    ) |>
    mutate(
        set_id = paste(
            cell_type,
            category |>
                str_replace_all('[()+-]', '') |>
                str_replace_all(' +', '_'),
            sep = '_'
        )
    ) |>
    arrange(set_id) |>
    select(set_id, gene_id) |>
    write_tsv(out_path)

session_info()
