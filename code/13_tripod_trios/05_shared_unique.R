library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)

cell_types = c(
    "Astrocyte", "Endo", "Ependymal", "Excit.Thal", "Inhib_LHb_4.1",
    "Inhib_LHb_4.2", "Inhib.Thal", "LHb.1.3.4", "LHb.2.7", "LHb.4", "MHb.1",
    "MHb.1.2", "MHb.2", "MHb.3", "Microglia", "Oligo", "OPC"
)
trio_paths = here(
    "processed-data", "13_tripod_trios", "03_tripod_trios", "trios_%s.parquet"
)
link_path = here(
    'processed-data', '12_new_peaks', '15_donor_specificity',
    'metacell_filtered_unbiased.parquet'
)
plot_dir = here("plots", "13_tripod_trios", "05_shared_unique")
FDR_thres = 0.05
cell_type_colors = c(
    "shared" = "gray",
    "Oligo" = "#4d5802",
    "OPC" = "#d3c871",
    "Microglia" = "#222222",
    "Astrocyte" = "#8d363c",
    "Endo" = "#ee6c14",
    "Inhib.Thal" = '#b5a2ff',
    "Excit.Thal" = "#9e4ad1",
    "LHb.1.3.4" = "#0085af",
    "LHb.2.7" = "#7DF9FF",
    "LHb.4" = '#76c2af',
    "MHb.1" = "#FF00FF",
    "MHb.1.2" = "#c76a6a",
    "MHb.2" = "#FAA0A0",
    "MHb.3" = "#fa246a",
    "Inhib_LHb_4.1" = "#23008b",
    "Inhib_LHb_4.2" = "#3333fa",
    "Ependymal" = "#f5a105ff"
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

order_cell_types = function(this_df) {
    this_df = this_df |>
        mutate(
            cell_type = factor(
                cell_type,
                levels = c(
                    this_df |>
                        group_by(cell_type) |>
                        summarise(n_links = n()) |>
                        arrange(desc(n_links)) |>
                        pull(cell_type) |>
                        as.character(),
                    setdiff(cell_types, as.character(unique(this_df$cell_type)))
                )
            )
        )
    return(this_df)
}

count_barplot = function(this_df, y_lab, plot_path) {
    p = ggplot(this_df, aes(x = cell_type, fill = shared_type)) +
        geom_bar() +
        scale_fill_manual(values = cell_type_colors) +
        scale_x_discrete(drop = FALSE) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        guides(fill = "none") +
        labs(x = "Cell Type", y = y_lab)
    pdf(plot_path, height = 5)
    print(p)
    dev.off()
}

################################################################################
#   Main
################################################################################

trio_df_list = list()
for (cell_type in cell_types) {
    this_path = sprintf(trio_paths, cell_type)
    if (file.exists(this_path)) {
        trio_df_list[[cell_type]] = read_parquet_duckdb(
                this_path, prudence = 'lavish'
            ) |>
            filter(
                condition_on == 'Yj', stringency_level == 2, adj < FDR_thres
            ) |>
            select(peak, gene, TF, adj) |>
            mutate(cell_type = cell_type)
    } else {
        message(
            sprintf(
                "Trios were not found for cell type %s; skipping", cell_type
            )
        )
    }
}
trio_df = bind_rows(trio_df_list) |>
    group_by(peak, gene, TF) |>
    mutate(
        shared_type = factor(
            ifelse(n() > 1, "shared", cell_type),
            levels = names(cell_type_colors)
        )
    ) |>
    ungroup() |>
    order_cell_types() |>
    collect()

link_df = read_parquet_duckdb(link_path) |>
    mutate(
        shared_type = factor(
            ifelse(is_shared, "shared", as.character(cell_type)),
            levels = names(cell_type_colors)
        )
    ) |>
    order_cell_types() |>
    collect()

shared_df = inner_join(
        trio_df |>
            select(-shared_type),
        link_df |>
            filter(score > 0),
        by = c("peak", "gene", "cell_type")
    ) |>
    order_cell_types()

count_barplot(trio_df, "Number of Trios", file.path(plot_dir, "trios.pdf"))
count_barplot(
    shared_df, "# of Trio-Link Overlaps",
    file.path(plot_dir, "overlaps.pdf")
)

session_info()
