library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    "filtered_trios_fine.parquet"
)
plot_dir = here("plots", "13_tripod_trios", "14_TF_composition")

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

trio_df = read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    filter(is_intersect, stringency_level == 1) |>
    select(TF, cell_type) |>
    collect()

# ── TF composition stacked bar plot ──────────────────────────────────────────

top20_tfs <- trio_df |>
    count(TF, sort = TRUE) |>
    slice_head(n = 20) |>
    pull(TF)

palette_20 <- c(
    "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4",
    "#42D4F4", "#F032E6", "#BFEF45", "#FABED4", "#469990",
    "#DCBEFF", "#9A6324", "#FFFAC8", "#800000", "#AAFFC3",
    "#808000", "#FFD8B1", "#000075", "#A9A9A9", "#000000"
)
tf_colors <- setNames(palette_20, top20_tfs)

p_composition <- trio_df |>
    mutate(
        TF_grouped = if_else(TF %in% top20_tfs, TF, "Other"),
        TF_grouped = factor(TF_grouped, levels = c(top20_tfs, "Other"))
    ) |>
    count(cell_type, TF_grouped) |>
    group_by(cell_type) |>
    mutate(prop = n / sum(n)) |>
    ungroup() |>
    ggplot(aes(x = cell_type, y = prop, fill = TF_grouped)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(tf_colors, "Other" = "gray80")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme_bw(base_size = 18) +
    labs(x = "Cell Type", y = "Proportion of Trios", fill = "TF")

pdf(file.path(plot_dir, "TF_composition_by_cell_type.pdf"))
print(p_composition)
dev.off()

session_info()
