library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(cowplot)

in_dir = here('processed-data', '15_DARs', '06_cell_level_DARs')
out_dir = here('processed-data', '15_DARs', '07_cell_level_gather')
plot_dir = here('plots', '15_DARs', '07_cell_level_gather')
cell_map_path = here("raw-data", "cell_type_map.csv")
colors_path = here(
    'code', '05_03_annotation_adjustments', 'celltype_colors.R'
)
p_adj_cutoff = 0.05
log_fc_cutoff = 1
max_dars = 10000

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

source(colors_path)
cell_type_colors = my_colors_mid

cluster_map = read_csv(cell_map_path, show_col_types = FALSE)
rename_map = stats::setNames(
    cluster_map$new_cell_type, cluster_map$old_cell_type
)

################################################################################
#   Functions
################################################################################

my_barplot = function(dar_df) {
    cell_types = unique(dar_df$cell_type)

    dar_df = dar_df |>
        group_by(cell_type) |>
        summarise(num_DARs = n(), .groups = 'drop') |>
        right_join(
            tibble(cell_type = cell_types),
            by = 'cell_type'
        ) |>
        mutate(
            num_DARs = replace_na(num_DARs, 0),
            cell_type = dplyr::coalesce(
                    unname(rename_map[as.character(cell_type)]),
                    as.character(cell_type)
                ) |>
                factor(levels = cluster_map$new_cell_type)
        )
    
    cell_types = dplyr::coalesce(
        unname(rename_map[as.character(cell_types)]),
        as.character(cell_types)
    )

    p = dar_df |>
        ggplot(aes(x = cell_type, y = num_DARs, fill = cell_type)) +
            geom_bar(stat = 'identity') +
            scale_fill_manual(values = cell_type_colors[cell_types]) +
            scale_y_log10() +
            theme_bw(base_size = 18) +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            ) +
            guides(fill = 'none') +
            labs(x = "Cell type", y = "Number of DARs")
    
    return(p)
}

################################################################################
#   Read in and perform some filtering
################################################################################

in_files = list.files(
    in_dir, pattern = "\\.parquet$", full.names = TRUE
)
dar_df_list = list()
for (in_file in in_files) {
    dar_df_list[[in_file]] = read_parquet_duckdb(
            in_file, prudence = "lavish"
        ) |>
        select(peak, cell_type, avg_log2FC, p_val_adj) |>
        filter(p_val_adj < p_adj_cutoff, abs(avg_log2FC) > log_fc_cutoff)
}
dar_df = bind_rows(dar_df_list) |>
    collect()

################################################################################
#   Plot and export
################################################################################

pdf(file.path(plot_dir, "DARs_per_cell_type.pdf"), width = 10, height = 5)
print(my_barplot(dar_df))
dev.off()

#   For LDSC, limit the number of DARs and retain just basic info
dar_df |>
    mutate(DA_direction = ifelse(sign(avg_log2FC) > 0, "Up", "Down")) |>
    group_by(cell_type, DA_direction) |>
    arrange(desc(abs(avg_log2FC))) |>
    slice_head(n = max_dars) |>
    ungroup() |>
    select(peak, cell_type, DA_direction) |>
    arrange(cell_type, DA_direction, peak) |>
    write_csv(file.path(out_dir, "DARs_LDSC.csv.gz"))

#   For supp table and Shiny app, include all significant DARs
write_csv(dar_df, file.path(out_dir, "DARs_all.csv.gz"))

session_info()
