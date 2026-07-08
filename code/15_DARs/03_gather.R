library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(cowplot)

in_dir = here('processed-data', '15_DARs', '02_calculate_DARs')
out_dir = here('processed-data', '15_DARs', '03_gather')
plot_dir = here('plots', '15_DARs', '03_gather')
task_map_path = here(
    'processed-data', '15_DARs', '01_pseudobulk_atac',
    'task_map.csv'
)
colors_path = here(
    'code', '05_03_annotation_adjustments', 'celltype_colors.R'
)
p_adj_cutoff = 0.05
log_fc_cutoff = log2(1.5)
max_dars = 1000

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

source(colors_path)
cell_type_colors = c(my_colors_mid, my_colors_class)
cell_type_colors[['Neuron']] = '#532222'

################################################################################
#   Functions
################################################################################

my_barplot = function(dar_df, cell_types, resolution, lab_title) {
    p = dar_df |>
        group_by(cell_type) |>
        summarise(num_DARs = n(), .groups = 'drop') |>
        right_join(
            tibble(cell_type = cell_types),
            by = 'cell_type'
        ) |>
        mutate(num_DARs = replace_na(num_DARs, 0)) |>
        mutate(cell_type = factor(cell_type, levels = cell_types)) |>
        ggplot(aes(x = cell_type, y = num_DARs, fill = cell_type)) +
            geom_bar(stat = 'identity') +
            scale_fill_manual(values = cell_type_colors[cell_types]) +
            scale_y_log10() +
            theme_bw() +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            ) +
            guides(fill = 'none') +
            labs(x = "Cell type", y = "Number of DARs", title = lab_title)
    
    return(p)
}

plot_and_process = function(dar_df, task_map_df, this_resolution) {
    these_cell_types = task_map_df |>
        filter(resolution == this_resolution) |>
        pull(cell_type)

    this_dar_df = dar_df |>
        filter(resolution == this_resolution)

    p_before = my_barplot(
        this_dar_df, these_cell_types, this_resolution, lab_title = "Before"
    )

    this_dar_df = this_dar_df |>
        filter(avg_log2FC > log_fc_cutoff) |>
        group_by(peak) |>
        filter(n() == 1) |>
        ungroup()

    p_after = my_barplot(
        this_dar_df, these_cell_types, this_resolution, lab_title = "After"
    )
    
    pdf(file.path(plot_dir, sprintf("DAR_barplot_%s.pdf", this_resolution)))
    print(plot_grid(p_before, p_after, ncol = 1))
    dev.off()

    this_dar_df |>
        group_by(cell_type) |>
        arrange(desc(avg_log2FC)) |>
        slice_head(n = max_dars) |>
        ungroup() |>
        write_csv(
            file.path(out_dir, sprintf("DARs_%s.csv.gz", this_resolution))
        )
}

################################################################################
#   Read in and perform minimal filtering
################################################################################

in_files = list.files(
    in_dir, pattern = "\\.parquet$", full.names = TRUE
)
dar_df_list = list()
for (in_file in in_files) {
    dar_df_list[[in_file]] = read_parquet_duckdb(
            in_file, prudence = "lavish"
        ) |>
        select(peak, cell_type, resolution, avg_log2FC, p_val_adj) |>
        filter(p_val_adj < p_adj_cutoff, avg_log2FC > 0)
}
dar_df = bind_rows(dar_df_list) |>
    collect()

################################################################################
#   Group into resolutions as we'll use downstream
################################################################################

dar_df_list = list()

dar_df_list[['broad']] = dar_df |>
    filter(!str_detect(cell_type, '[ML]Hb|Thal')) |>
    mutate(resolution = 'broad')

dar_df_list[['mid']] = dar_df |>
    filter(
        (resolution == 'mid') |
        ((resolution == 'fine') & !str_detect(cell_type, '[ML]Hb|Thal'))
    ) |>
    mutate(resolution = 'mid')

dar_df_list[['fine']] = dar_df |>
    filter(resolution == 'fine') |>
    mutate(resolution = 'fine')

dar_df = bind_rows(dar_df_list)

task_map = read_csv(task_map_path, show_col_types = FALSE)
task_map_df_list = list()

task_map_df_list[['broad']] = task_map |>
    filter(!str_detect(cell_type, '[ML]Hb|Thal')) |>
    mutate(resolution = 'broad')

task_map_df_list[['mid']] = task_map |>
    filter(
        (resolution == 'mid') |
        ((resolution == 'fine') & !str_detect(cell_type, '[ML]Hb|Thal'))
    ) |>
    mutate(resolution = 'mid')

task_map_df_list[['fine']] = task_map |>
    filter(resolution == 'fine') |>
    mutate(resolution = 'fine')

task_map_df = bind_rows(task_map_df_list) |>
    select(resolution, cell_type)

################################################################################
#   Plot and export for each resolution
################################################################################

for (this_resolution in c('broad', 'mid', 'fine')) {
    plot_and_process(dar_df, task_map_df, this_resolution)
}

session_info()
