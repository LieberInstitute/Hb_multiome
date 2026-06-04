library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
plot_dir = here("plots", "13_tripod_trios", "12_TF_overlap")
cell_types = c("MHb.2", "LHb.2.7")

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

trio_df = read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    select(peak, gene, TF, cell_type, stringency_level) |>
    collect()

top_trio_df = read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    filter(cell_type %in% cell_types, is_intersect, is_unique, is_top_TF) |>
    collect() |>
    group_by(cell_type, peak, gene, TF) |>
    filter(n() == 2) |>
    ungroup() |>
    distinct(cell_type, peak, gene, TF)

#-------------------------------------------------------------------------------
#   TF breadth: histogram of how many cell types each TF appears in
#-------------------------------------------------------------------------------

p = trio_df |>
    distinct(TF, cell_type, stringency_level) |>
    count(TF, stringency_level, name = "n_cell_types") |>
    ggplot(aes(x = n_cell_types)) +
        geom_histogram(binwidth = 1) +
        facet_wrap(~stringency_level, nrow = 2, labeller = label_both) +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        scale_y_continuous(labels = scales::comma) +
        labs(x = "Number of cell types", y = "Number of TFs") +
        theme_bw(base_size = 15)
pdf(file.path(plot_dir, "TF_breadth_histogram.pdf"))
print(p)
dev.off()

#-------------------------------------------------------------------------------
#   Top 5 TFs per cell type by number of trios
#-------------------------------------------------------------------------------

for (lvl in 1:2) {
    p = trio_df |>
        filter(stringency_level == lvl) |>
        count(cell_type, TF, name = "n_trios") |>
        group_by(cell_type) |>
        slice_max(n_trios, n = 5, with_ties = FALSE) |>
        ungroup() |>
        mutate(
            TF = paste(cell_type, TF, sep = "___"),
            TF = reorder(TF, n_trios)
        ) |>
        ggplot(aes(x = n_trios, y = TF)) +
        geom_col() +
        scale_y_discrete(labels = \(x) sub(".*___", "", x)) +
        facet_wrap(~cell_type, scales = "free") +
        labs(
            x = "Number of trios", y = NULL,
            title = sprintf("Top 5 TFs per cell type (stringency level %d)", lvl)
        ) +
        theme_bw(base_size = 13)

    pdf(file.path(plot_dir, sprintf("Top_5_TFs_level_%d.pdf", lvl)))
    print(p)
    dev.off()
}

#-------------------------------------------------------------------------------
#   Number of cell types each TF in top_trio_df appears in across all trios
#-------------------------------------------------------------------------------

message("Cell-type specificty of TFs from top trios:")
trio_df |>
    filter(TF %in% top_trio_df$TF) |>
    distinct(TF, cell_type) |>
    count(TF, name = "n_cell_types") |>
    arrange(desc(n_cell_types)) |>
    print()

session_info()
