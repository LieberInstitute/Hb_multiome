library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)

cell_type_res = c("broad", "fine")[
    as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
]
trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    sprintf("filtered_trios_%s.parquet", cell_type_res)
)
plot_dir = here("plots", "13_tripod_trios", "05_shared_unique")

if (cell_type_res == 'broad') {
    my_colors_mid = c(
        shared = "gray",
        MHb = "#ad1d8c",
        LHb = "#1f78b4",
        Inhib_LHb = "#c70404"
    )
} else {
    source(here("code", "05_03_annotation_adjustments", "celltype_colors.R"))
    my_colors_mid[['shared']] = "gray"
}

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
        scale_fill_manual(values = my_colors_mid) +
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

trio_df = read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    collect() |>
    mutate(
        shared_type = factor(
            ifelse(is_unique, cell_type, "shared"),
            levels = names(my_colors_mid)
        )
    ) |>
    order_cell_types()

for (this_stringency_level in c(1, 2)) {
    count_barplot(
        trio_df |> filter(stringency_level == this_stringency_level),
        "Number of Trios",
        file.path(
            plot_dir,
            sprintf(
                "trios_level%d_%s.pdf", this_stringency_level, cell_type_res
            )
        )
    )
}

session_info()
