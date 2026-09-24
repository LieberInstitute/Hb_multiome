#   Do SCZD DEGs from habenula pilot overlap any genes linked to peaks in this
#   project?

library(tidyverse)
library(here)
library(sessioninfo)

deg_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/10_DEA/04_DEA/DEA_All-gene_qc-totAGene-qSVs-Hb-Thal.tsv'
peak_path = here(
    'processed-data', '06_peak_calling', '21_overlaping_FDRscores_TopHeatmap',
    'all_peaks_categorized.csv.gz'
)

peak_df = read_csv(peak_path, show_col_types = FALSE) |>
    filter(!is.na(link_gene_id)) |>
    select(cell_type, link_gene_id, link_gene_name, category)

deg_df = read_tsv(deg_path, show_col_types = FALSE) |>
    filter(adj.P.Val < 0.1) |>
    dplyr::rename(
        link_gene_id = ensemblID, SCZD_FDR = adj.P.Val, SCZD_logFC = logFC
    ) |>
    select(link_gene_id, SCZD_FDR, SCZD_logFC)

shared_df = deg_df |>
    inner_join(peak_df, by = 'link_gene_id')

message("All SCZD DEGs linked to peaks at FDR < 0.1:")
print(shared_df, n = nrow(shared_df))
    
message("All SCZD DEGs linked to peaks at FDR < 0.05:")
temp = shared_df |> filter(SCZD_FDR < 0.05)
print(temp, n = nrow(temp))

session_info()
