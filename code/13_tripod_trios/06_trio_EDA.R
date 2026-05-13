library(sessioninfo)
library(tidyverse)
library(here)
library(duckplyr)
library(TRIPOD)
library(qs2)
library(cowplot)

source(here("code", "05_03_annotation_adjustments", "celltype_colors.R"))

cell_types = c(
    "Astrocyte", "Endo", "Ependymal", "Excit.Thal", "Inhib_LHb_4.1",
    "Inhib_LHb_4.2", "Inhib.Thal", "LHb.1.3.4", "LHb.2.7", "LHb.4", "MHb.1",
    "MHb.1.2", "MHb.2", "MHb.3", "Microglia", "Oligo", "OPC"
)
trio_paths = here(
    "processed-data", "13_tripod_trios", "03_tripod_trios", "trios_%s.parquet"
)
link_path = here(
    'processed-data', '12_new_peaks', '15_donor_specificity',
    'metacell_filtered_unbiased.parquet'
)
model_paths = here(
    "processed-data", "13_tripod_trios", "03_tripod_trios", "fit_models_%s.qs2"
)
plot_dir = here("plots", "13_tripod_trios", "06_trio_EDA")
num_examples = 5

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Functions
################################################################################

plot_facet_by_set = function(set_1, set_2, y_lab, plot_path) {
    intersect_set = inner_join(
        set_1, set_2, by = setdiff(names(set_1), "category")
    )
  
    p = bind_rows(
            set_1,
            set_2,
            intersect_set |> mutate(category = "Intersection")
        ) |>
        count(cell_type, category) |>
        mutate(
            cell_type = factor(cell_type, levels = cell_types),
            category = factor(category, levels = c(
                unique(set_1$category), unique(set_2$category), "Intersection"
            ))
        ) |>
        ggplot(aes(x = cell_type, y = n, fill = cell_type)) +
            geom_col() +
            scale_y_continuous(labels = scales::comma) +
            scale_fill_manual(values = my_colors_mid) +
            facet_wrap(~category, nrow = 3, scales = "free_y") +
            labs(x = NULL, y = y_lab) +
            guides(fill = "none") +
            theme_bw(base_size = 15) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    pdf(plot_path)
    print(p)
    dev.off()

    invisible(NULL)
}

plot_facet_by_celltype = function(set_1, set_2, y_lab, plot_path) {
    intersect_set = inner_join(
        set_1, set_2, by = setdiff(names(set_1), "category")
    )
  
    p = bind_rows(
            set_1,
            set_2,
            intersect_set |> mutate(category = "Intersection")
        ) |>
        count(cell_type, category) |>
        mutate(category = factor(category, levels = c(
            unique(set_1$category), unique(set_2$category), "Intersection"
        ))) |>
        ggplot(aes(x = category, y = n, fill = cell_type)) +
            geom_col() +
            facet_wrap(~cell_type, scales = "fixed") +
            scale_y_continuous(labels = scales::comma, transform = "log10") +
            scale_fill_manual(values = my_colors_mid) +
            labs(x = NULL, y = y_lab) +
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

trio_df_list = list()
for (cell_type in cell_types) {
    this_path = sprintf(trio_paths, cell_type)
    if (file.exists(this_path)) {
        trio_df_list[[cell_type]] = read_parquet_duckdb(
                this_path, prudence = 'lavish'
            ) |>
            filter(condition_on == 'Yj', adj < 0.05) |>
            select(peak, gene, TF, coef, adj, stringency_level) |>
            mutate(cell_type = cell_type)
    } else {
        message(
            sprintf(
                "Trios were not found for cell type %s; skipping", cell_type
            )
        )
    }
}
trio_df = bind_rows(trio_df_list) |>
    collect()

link_df = read_parquet_duckdb(link_path, prudence = 'stingy') |>
    select(peak, gene, cell_type) |>
    collect()

#-------------------------------------------------------------------------------
#   Trio-link overlap
#-------------------------------------------------------------------------------

trio_set = trio_df |>
    filter(stringency_level == 1) |>
    distinct(peak, gene, cell_type) |>
    mutate(category = "Trios")

link_set = link_df |>
    distinct(peak, gene, cell_type) |>
    mutate(category = "Linked peaks")

plot_facet_by_set(
    trio_set, link_set,
    y_lab   = "Number of peak-gene pairs",
    plot_path = file.path(plot_dir, "trio_link_overlap_by_set.pdf")
)

plot_facet_by_celltype(
    trio_set, link_set,
    y_lab   = "Number of peak-gene pairs",
    plot_path = file.path(plot_dir, "trio_link_overlap_by_celltype.pdf")
)

#-------------------------------------------------------------------------------
#   Trio-stringency overlap
#-------------------------------------------------------------------------------

level1_set = trio_df |>
    filter(stringency_level == 1) |>
    distinct(peak, gene, TF, cell_type) |>
    mutate(category = "Level 1")

level2_set = trio_df |>
    filter(stringency_level == 2) |>
    distinct(peak, gene, TF, cell_type) |>
    mutate(category = "Level 2")

plot_facet_by_set(
    level1_set, level2_set,
    y_lab   = "Number of significant trios",
    plot_path = file.path(plot_dir, "trio_counts_by_set.pdf")
)

plot_facet_by_celltype(
    level1_set, level2_set,
    y_lab   = "Number of significant trios",
    plot_path = file.path(plot_dir, "trio_counts_by_celltype.pdf")
)

xy_mat_list = lapply(
    unique(trio_df$cell_type),
    function(this_cell_type) {
        model_paths |>
            sprintf(this_cell_type) |>
            qs_read()
    }
)
names(xy_mat_list) = unique(trio_df$cell_type)

for (this_test_level in 1:2) {
    plot_list = list()
    for (this_cell_type in names(xy_mat_list)) {
        sub_df = trio_df |>
            filter(
                cell_type == this_cell_type,
                stringency_level == this_test_level
            ) |>
            slice_sample(n = num_examples)

        for (i in seq_len(num_examples)) {
            plot_list[[length(plot_list) + 1]] = plotGenePeakTFScatter(
                xymats = xy_mat_list[[this_cell_type]][[sub_df$gene[i]]],
                peak.name = sub_df$peak[i], TF.name = sub_df$TF[i],
                to.plot = "TRIPOD", match.by = "Yj",
                level = this_test_level, cap.at.quantile = 0
            ) + labs(title = sprintf("Example %d", i)) #+
            #   theme_bw(base_size = 15)
        }
    }

    pdf(file.path(plot_dir, sprintf("trio_scatter_level_%d.pdf", this_test_level)))
    print(plot_grid(plotlist = plot_list, ncol = length(xy_mat_list)))
    dev.off()
}

session_info()
