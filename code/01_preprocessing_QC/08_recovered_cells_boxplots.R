library(here)
library(tidyverse)
library(sessioninfo)

name_map = tibble(
    run_name = c("ARC", "GEX", "ATAC", "ARC_reanalyze"),
    run_label = c(
        "Cellranger-ARC",
        "Cellranger-count",
        "Cellranger-ATAC",
        "Cellranger-ARC reanalyze"
    )
)
plot_dir = here('plots', '01_preprocessing_QC', '08_recovered_cells_boxplots')

dir.create(plot_dir, showWarnings = FALSE)

plot_df_list = list()
for (run_name in c("GEX", "ARC", "ATAC")) {
    csv_files = list.files(
        here(
            'processed-data' , sprintf('cellranger%s_summary_rpts', run_name),
            'csv_web_summary_rpt'
        ),
        pattern = '^S([3-9]|[0-9]{2}).*_summary\\.csv$',
        full.names = TRUE
    )

    plot_df_list[[run_name]] = tibble(
        num_cells = sapply(
            csv_files,
            function(x) {
                read_csv(x, show_col_types = FALSE) |>
                    pull(matches("[Ee]stimated [Nn]umber of [Cc]ells"))
            }
        ),
        run_name = run_name
    )
}

csv_files = list.files(
    here(
        'processed-data', '01_preprocessing_QC', 'cellrangerARC_reanalyze',
        'csv_files'
    ),
    pattern = '_total_TrueCells\\.csv$',
    full.names = TRUE
)
stopifnot(length(csv_files) == 10)

plot_df_list[['ARC_reanalyze']] = tibble(
    num_cells = sapply(
        csv_files,
        function(x) {
            read_csv(x, show_col_types = FALSE) |>
                filter(Cells == 'NonEmptyCells') |>
                pull(values) |>
                as.numeric()
        }
    ) |> unname(),
    run_name = 'ARC_reanalyze'
)

p = bind_rows(plot_df_list) |>
    left_join(name_map, by = 'run_name') |>
    mutate(run_label = factor(run_label, levels = name_map$run_label)) |>
    ggplot(aes(x = run_label, y = num_cells)) +
        geom_boxplot() +
        labs(x = "Pipeline", y = "Recovered Cells") +
        theme_bw(base_size = 18) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
pdf(file.path(plot_dir, 'recovered_cells.pdf'), width = 5)
print(p)
dev.off()

session_info()
