#   Using the data from every combination of two cell types, explore how the
#   distribution of correlation scores differ between the target and non-target
#   cell types

library(here)
library(tidyverse)
library(sessioninfo)
library(duckplyr)

cell_types = c(
    "Astrocyte", "Endo", "Ependymal", "Excit.Thal", "Inhib_LHb_4.1",
    "Inhib_LHb_4.2", "Inhib.Thal", "LHb.1.3.4", "LHb.2.7", "LHb.4", "MHb.1",
    "MHb.1.2", "MHb.2", "MHb.3", "Microglia", "Oligo", "OPC"
)
link_path = here(
    'processed-data', '12_new_peaks', '04_filter_links',
    'metacell_all_data.parquet'
)
plot_dir = here("plots", "12_new_peaks", "10_off_target_cor")
fdr_cutoffs = c(0.01, 0.05, 0.1, 0.2, 1)
num_rows_sample = 1e6

set.seed(0)
num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

result_df = read_parquet_duckdb(link_path, prudence = 'stingy') |>
    select(peak_called, score, FDR) |>
    collect()

#   Density plots of correlation scores split by whether target cell type matches
#   other cell type. No substantial difference is seen at this level
p = ggplot(
        result_df,
        aes(x = abs(score), fill = peak_called, color = peak_called)
    ) +
    geom_density(alpha = 0.5, linewidth = 1) +
    theme_bw(base_size = 20) +
    labs(
        x = "abs(Correlation Score)", y = "Density", fill = "Cell Type",
        color = "Cell Type"
    )
pdf(file.path(plot_dir, "off_target_cor_global.pdf"), width = 10, height = 5)
print(p)
dev.off()

result_df_list = list()
for (fdr_cutoff in fdr_cutoffs) {
    result_df_list[[as.character(fdr_cutoff)]] = result_df |>
        filter(FDR <= fdr_cutoff) |>
        mutate(FDR_class = sprintf("FDR <= %s", fdr_cutoff)) |>
        #   There are too many rows (> 9e8). Sample to control memory
        slice_sample(n = num_rows_sample)
}
result_df = bind_rows(result_df_list) |>
    mutate(
        FDR_class = factor(
            FDR_class, levels = sprintf("FDR <= %s", fdr_cutoffs)
        )
    )

#   Stratifying by FDR class also shows no substantial difference in the
#   distributions
p = ggplot(
        result_df,
        aes(x = abs(score), fill = peak_called, color = peak_called)
    ) +
    geom_density(alpha = 0.5, linewidth = 1) +
    facet_wrap(~ FDR_class, nrow = 4) +
    labs(
        x = "abs(Correlation Score)", y = "Density", fill = "Target Cell Type",
        color = "Target Cell Type"
    ) +
    theme_bw(base_size = 20)
pdf(file.path(plot_dir, "off_target_cor_stratified.pdf"), width = 10, height = 10)
print(p)
dev.off()

session_info()
