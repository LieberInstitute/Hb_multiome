library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(TRIPOD)
library(qs2)
library(cowplot)

cell_type_res = c('broad', 'fine')[
    as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
]

if (cell_type_res == 'broad') {
    cell_types = c(
        "Astrocyte", "Endo", "Ependymal", "Excit.Thal", "Inhib.Thal",
        "Microglia", "Oligo", "OPC", "MHb", "LHb", "Inhib_LHb"
    )

    my_colors_mid = c(
        shared = "gray",
        MHb = "#ad1d8c",
        LHb = "#1f78b4",
        Inhib_LHb = "#c70404"
    )
} else {
    cell_types = c(
        "Astrocyte", "Endo", "Ependymal", "Excit.Thal", "Inhib_LHb_4.1",
        "Inhib_LHb_4.2", "Inhib.Thal", "LHb.1.3.4", "LHb.2.7", "LHb.4", "MHb.1",
        "MHb.1.2", "MHb.2", "MHb.3", "Microglia", "Oligo", "OPC"
    )
    cell_type_levels = c(
        "MHb_A", "MHb_B", "MHb_C", "MHb_D", "LHb_A", "LHb_B", "LHb_C",
        "GABA_LHb_C.1", "GABA_LHb_C.2", "Excit.Thal", "Inhib.Thal",
        "Astrocyte", "Endo", "Ependymal", "Microglia", "Oligo", "OPC"
    )
    cell_map_path = here("raw-data", "cell_type_map.csv")

    source(here("code", "05_03_annotation_adjustments", "celltype_colors.R"))
    my_colors_mid[['shared']] = "gray"
    my_colors_mid = my_colors_mid[c('shared', cell_type_levels)]
}

trio_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    sprintf("filtered_trios_%s.parquet", cell_type_res)
)
model_paths = here(
    "processed-data", "13_tripod_trios", "03_tripod_trios", "fit_models_%s.qs2"
)
out_path = here(
    "processed-data", "13_tripod_trios", "06_trio_EDA",
    sprintf("trio_summary_%s.csv", cell_type_res)
)
plot_dir = here("plots", "13_tripod_trios", "06_trio_EDA")
num_examples = 5
num_trios_export = 50

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)
dir.create(dirname(out_path), showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

plot_cell_type = identity
if (cell_type_res == "fine") {
    cluster_map = read_csv(cell_map_path, show_col_types = FALSE)
    rename_map = stats::setNames(cluster_map$new_cell_type, cluster_map$old_cell_type)
    plot_cell_type = function(x) dplyr::coalesce(unname(rename_map[x]), x)
}

plot_facet_by_set = function(intersect_set, plot_path) {
    p = ggplot(intersect_set, aes(x = cell_type_plot, fill = shared_type)) +
        geom_bar() +
        scale_y_continuous(labels = scales::comma) +
        scale_fill_manual(values = my_colors_mid) +
        facet_wrap(~category, nrow = 3, scales = "free_y") +
        labs(x = NULL, y = "Number of significant trios") +
        guides(fill = "none") +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    pdf(plot_path)
    print(p)
    dev.off()

    invisible(NULL)
}

plot_facet_by_celltype = function(intersect_set, plot_path) {
    p = ggplot(intersect_set, aes(x = category, fill = cell_type_plot)) +
        geom_bar() +
        facet_wrap(~cell_type_plot, scales = "fixed") +
        scale_y_continuous(transform = "log10", labels = scales::comma) +
        scale_fill_manual(values = my_colors_mid) +
        labs(x = NULL, y = "Number of significant trios") +
        guides(fill = "none") +
        theme_bw(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    pdf(plot_path)
    print(p)
    dev.off()
  
    invisible(NULL)
}

################################################################################
#   Main
################################################################################

trio_df = read_parquet_duckdb(trio_path, prudence = 'stingy') |>
    collect() |>
    mutate(
        cell_type_plot = plot_cell_type(cell_type),
        shared_type = factor(
            ifelse(is_unique_intersect | !is_intersect, cell_type_plot, "shared"),
            levels = names(my_colors_mid)
        ),
        cell_type_plot = factor(cell_type_plot, levels = setdiff(names(my_colors_mid), "shared")),
        cell_type = factor(cell_type, levels = cell_types)
    )

#-------------------------------------------------------------------------------
#   Export easy-to-view CSV of top trios for each cell type and test level
#-------------------------------------------------------------------------------

overlap_df = trio_df |>
    filter(stringency_level == 1, is_intersect) |>
    mutate(stringency_level = "Intersection") |>
    select(peak, gene, TF, cell_type, adj, stringency_level, is_unique_intersect)

trio_df |>
    select(peak, gene, TF, cell_type, adj, stringency_level, is_unique_intersect) |>
    rbind(overlap_df) |>
    dplyr::rename(p_adj = adj, TRIPOD_test_level = stringency_level) |>
    select(
        cell_type, TRIPOD_test_level, peak, gene, TF, p_adj, is_unique_intersect
    ) |>
    group_by(cell_type, TRIPOD_test_level) |>
    arrange(p_adj) |>
    slice_head(n = num_trios_export) |>
    ungroup() |>
    arrange(TRIPOD_test_level, cell_type, is_unique_intersect, p_adj) |>
    write_csv(out_path)

#-------------------------------------------------------------------------------
#   Level 1 and level 2 overlap
#-------------------------------------------------------------------------------

level1_set = trio_df |>
    filter(stringency_level == 1) |>
    select(peak, gene, TF, cell_type, cell_type_plot, shared_type) |>
    mutate(category = "Level 1")

level2_set = trio_df |>
    filter(stringency_level == 2) |>
    select(peak, gene, TF, cell_type, cell_type_plot, shared_type) |>
    mutate(category = "Level 2")

intersect_set = inner_join(
        level1_set |> select(-category),
        level2_set |> select(-c(category, shared_type)),
        by = c("peak", "gene", "TF", "cell_type", "cell_type_plot")
    ) |>
    mutate(category = "Intersection")
  
all_data = bind_rows(level1_set, level2_set, intersect_set) |>
    mutate(
        category = factor(
            category, levels = c("Level 1", "Level 2", "Intersection")
        )
    )

plot_facet_by_set(
    all_data, plot_path = file.path(plot_dir, "trio_counts_by_set.pdf")
)

plot_facet_by_celltype(
    all_data,
    plot_path = file.path(plot_dir, "trio_counts_by_celltype.pdf")
)

p = trio_df |>
    filter(is_intersect, stringency_level == 1, grepl('Hb', cell_type)) |>
    ggplot(aes(x = cell_type_plot, fill = cell_type_plot)) +
        geom_bar() +
        scale_y_continuous(labels = scales::comma) +
        scale_fill_manual(values = my_colors_mid) +
        labs(x = "Cell Type", y = "Significant trios") +
        guides(fill = "none") +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pdf(file.path(plot_dir, "main_figure_subset.pdf"), height = 5)
print(p)
dev.off()

session_info()
