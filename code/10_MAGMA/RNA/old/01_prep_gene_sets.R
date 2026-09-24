library(tidyverse)
library(here)
library(sessioninfo)

for (i in c("broad","fine")) {
model_path = here(
    'processed-data', '10_MAGMA', 'RNA', 'registration_banksy',
    'modeling_results', paste0(i,'.rds')
)
out_path = here(
    'processed-data', '10_MAGMA', 'RNA', 'gene_sets',
    paste0('enrichment_markers_',i,'.tsv')
)

sig_cutoff = 0.05

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

readRDS(model_path)$enrichment |>
    as_tibble() |>
    select(ensembl, matches('^fdr_')) |>
    pivot_longer(
        cols = starts_with("fdr_"),
        names_to = "banksy",
        values_to = "fdr"
    ) |>
    filter(fdr < sig_cutoff) |>
    mutate(set_id = sub('fdr_', '', banksy)) |>
    dplyr::rename(gene_id = ensembl) |>
    select(set_id, gene_id) |> 
    write_tsv(out_path)
}

session_info()