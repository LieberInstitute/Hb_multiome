library(here)
library(tidyverse)
library(sessioninfo)

csv_files = list.files(
    here(
        'processed-data' , 'cellrangerARC_summary_rpts',
        'csv_web_summary_rpt'
    ),
    pattern = '^S([3-9]|[0-9]{2}).*_summary\\.csv$',
    full.names = TRUE
)
out_path = here(
    'processed-data', '18_misc_supp_tables', '02_cellranger_arc_metrics',
    'cellranger_arc_metrics.csv'
)
wet_bench_df = tibble(
    donor = c(
        'Br8518', 'Br8651', 'Br8274', 'Br6522', 'Br6432', 'Br8667', 'Br8582',
        'Br8433', 'Br9017', 'Br9090'
    ),
    num_neuronal_nuclei_sorted = c(
        91783, 123550, 198325, 111383, 210021, 119292, 114842, 76578, 124028,
        192043
    ),
    num_glial_nuclei_sorted = c(
        648000, 200000, 200000, 500000, 500000, 500000, 321842, 489595,
        500000, 500000
    ),
    num_neurons_loaded = c(
        12950, 16987.5, 18787.5, 18225, 19068.8, 17718.8, 17062.5, 11025,
        18187.5, 18225
    ),
    num_glia_loaded = c(5800, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    concentration_post_nuclei_permeabilization_nuclei_per_uL = c(
        9375, 5662.5, 12525, 9112.5, 12712.5, 7087.5, 4875, 1837.5, 7275, 12150
    ),
    total_nuclei_loaded = c(
        18750, 16987.5, 18787.5, 18225, 19068.8, 17718.8, 17062.5, 11025,
        18187.5, 18225
    )
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

metric_df = lapply(csv_files, read_csv) |>
    bind_rows()

session_info()
