#   How many unique gene-peak pairs do we get for various FDR and correlation
#   thresholds?

library(tidyverse)
library(here)
library(sessioninfo)
library(duckplyr)

link_path = here(
    'processed-data', '12_new_peaks', '04_filter_links',
    'metacell_all_data.parquet'
)
plot_path = here(
    "plots", "12_new_peaks", "09_threshold_heatmap", "heatmap.pdf"
)
fdr_thresholds = c(0.01, 0.05, 0.1, 0.2, 1)
cor_thresholds = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(plot_path), showWarnings = FALSE)

link_df = read_parquet_duckdb(link_path, prudence = 'stingy') |>
    collect()

## Count unique peak-gene pairs passing each threshold combination
count_df = expand.grid(fdr = fdr_thresholds, cor = cor_thresholds) |>
    as_tibble() |>
    mutate(
        n_pairs = map2_int(fdr, cor, \(f, c) {
            link_df |>
                filter(FDR < f, abs(score) > c) |>
                distinct(peak, gene) |>
                nrow()
        })
    )

## Build heatmap
p = count_df |>
    mutate(
        fdr = factor(fdr),
        cor = factor(cor)
    ) |>
    ggplot(aes(x = fdr, y = cor, fill = n_pairs)) +
    geom_tile() +
    geom_text(aes(label = scales::comma(n_pairs)), color = "white", size = 3) +
    scale_fill_viridis_c(
        trans = "log10",
        labels = scales::comma,
        name = "Unique\npeak-gene pairs"
    ) +
    labs(
        x = "FDR threshold",
        y = "Correlation threshold (|score|)",
        title = "Unique peak-gene pairs by threshold"
    ) +
    theme_bw(base_size = 15)

pdf(plot_path, height = 6)
print(p)
dev.off()

session_info()
