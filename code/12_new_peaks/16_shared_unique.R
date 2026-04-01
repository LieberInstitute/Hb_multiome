#   Create a barplot showing, for each cell type, how many significant links
#   are shared vs unique

library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)

link_path = here(
    'processed-data', '12_new_peaks', '01_link_peaks', 'all_data.parquet'
)
plot_dir = here("plots", "12_new_peaks", "16_shared_unique")
cor_thres = 0.3
FDR_thres = 0.1

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
    "MHb.3" = "#fa246a"
)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

link_df = read_parquet_duckdb(link_path) |>
    filter(
        score > abs(cor_thres), FDR < FDR_thres,
        target_cell_type == other_cell_type
    ) |>
    select(peak, gene, target_cell_type, FDR) |>
    #   Here I use a more relaxed definition of cell-type specificity. This is
    #   really one of two ways of doing this with the results that we have
    group_by(peak, gene) |>
    mutate(link_type = ifelse(n() == 1, target_cell_type, 'shared')) |>
    ungroup() |>
    select(target_cell_type, link_type) |>
    collect()

p = link_df |>
    mutate(
        target_cell_type = factor(
            target_cell_type,
            levels = link_df |>
                group_by(target_cell_type) |>
                summarise(n_links = n()) |>
                arrange(desc(n_links)) |>
                pull(target_cell_type)
        ),
        link_type = factor(link_type, levels = names(cell_type_colors))
    ) |>
    ggplot(aes(x = target_cell_type, fill = link_type)) +
        geom_bar() +
        scale_fill_manual(values = cell_type_colors) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        guides(fill = "none") +
        labs(x = "Cell Type", y = "Number of Links")
pdf(file.path(plot_dir, "shared_unique_links.pdf"), height = 5)
print(p)
dev.off()

session_info()
