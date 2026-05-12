#   Create a barplot showing, for each cell type, how many significant links
#   are shared vs unique

library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)

link_path = here(
    'processed-data', '12_new_peaks', '05_donor_specificity',
    'metacell_filtered_unbiased.parquet'
)
plot_dir = here("plots", "12_new_peaks", "06_shared_unique")

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
    Inhib_LHb_4.1 = "#23008b",
    Inhib_LHb_4.2 = "#3333fa",
    Ependymal = "#f5a105ff"
)

cell_types = names(cell_type_colors)[names(cell_type_colors) != "shared"]

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

link_df = read_parquet_duckdb(link_path) |>
    collect()

link_df = link_df |>
    mutate(
        cell_type = factor(
            cell_type,
            levels = c(
                link_df |>
                    group_by(cell_type) |>
                    summarise(n_links = n()) |>
                    arrange(desc(n_links)) |>
                    pull(cell_type),
                setdiff(cell_types, as.character(unique(link_df$cell_type)))
            )
        ),
        link_type = factor(
            ifelse(is_shared, "shared", as.character(cell_type)),
            levels = names(cell_type_colors)
        )
    )

p = ggplot(link_df, aes(x = cell_type, fill = link_type)) +
        geom_bar() +
        scale_fill_manual(values = cell_type_colors) +
        scale_x_discrete(drop = FALSE) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        guides(fill = "none") +
        labs(x = "Cell Type", y = "Number of Links")
pdf(file.path(plot_dir, "shared_unique_links.pdf"), height = 5)
print(p)
dev.off()

session_info()
