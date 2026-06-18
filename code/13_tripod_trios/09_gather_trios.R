library(sessioninfo)
library(tidyverse)
library(here)
library(qs2)
library(Seurat)
library(Signac)
library(duckplyr)
library(cowplot)

cell_type_res = c('broad', 'fine')[
    as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
]

if (cell_type_res == 'broad') {
    cell_types = c(
        "Astrocyte", "Endo", "Ependymal", "Excit.Thal", "Inhib.Thal",
        "Microglia", "Oligo", "OPC", "MHb", "LHb", "Inhib_LHb"
    )
    cell_types_keep = c("MHb", "LHb", "Inhib_LHb")
} else {
    cell_types = c(
        "Astrocyte", "Endo", "Ependymal", "Excit.Thal", "Inhib_LHb_4.1",
        "Inhib_LHb_4.2", "Inhib.Thal", "LHb.1.3.4", "LHb.2.7", "LHb.4", "MHb.1",
        "MHb.1.2", "MHb.2", "MHb.3", "Microglia", "Oligo", "OPC"
    )
    cell_types_keep = cell_types
}
trio_paths = here(
    "processed-data", "13_tripod_trios", "03_tripod_trios", "trios_%s.parquet"
)
prep_path = here(
    "processed-data", "13_tripod_trios", "02_tripod_preprocess",
    "preprocessed_objects_%s.qs2"
)
out_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    sprintf("filtered_trios_%s.parquet", cell_type_res)
)
scatter_data_path = here(
    "processed-data", "13_tripod_trios", "09_gather_trios",
    sprintf("scatter_plot_data_%s.rds", cell_type_res)
)
plot_dir = here("plots", "13_tripod_trios", "09_gather_trios", cell_type_res)
fdr_thres = 0.05
intersect_fdr_thres = 0.1
min_nonzero = 10
num_examples = 5

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(dirname(out_path), showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(0)

################################################################################
#   Functions
################################################################################

annotate_trios = function(df, bg_pgtf_list) {
    # bg_pgtf_list is a named list: names are cell types, values are unique
    # peak|gene|TF strings from sparsity-filtered FDR < intersect_fdr_thres trios
    # is_unique_intersect: TRUE only for intersect trios whose peak-gene-TF is
    # not present in any other cell type at either stringency level; FALSE
    # for all non-intersect trios
    df |>
        mutate(pg_tf = paste(peak, gene, TF, sep = "|")) |>
        group_by(cell_type) |>
        mutate(
            is_intersect = pg_tf %in% pg_tf[stringency_level == 1] &
                pg_tf %in% pg_tf[stringency_level == 2],
            is_unique_intersect = is_intersect & !pg_tf %in% unlist(
                bg_pgtf_list[names(bg_pgtf_list) != cell_type[1]]
            )
        ) |>
        group_by(peak, gene, cell_type, stringency_level) |>
        arrange(adj) |>
        mutate(is_top_TF = row_number() == 1) |>
        ungroup() |>
        arrange(cell_type, stringency_level, adj) |>
        select(-pg_tf)
}

make_scatter_data = function(trios_to_plot) {
    plot_data_list = list()
    for (ct in unique(trios_to_plot$cell_type)) {
        mats = qs_read(sprintf(prep_path, ct))$metacell_seur
        ct_trios = trios_to_plot |> filter(cell_type == ct)

        for (i in seq_len(nrow(ct_trios))) {
            row     = ct_trios[i, ]
            expr    = mats$rna[, row$gene, drop = TRUE]
            access  = mats$peak[, row$peak, drop = TRUE]
            tf_expr = mats$rna[, row$TF, drop = TRUE]

            plot_data_list[[length(plot_data_list) + 1]] = tibble(
                cell_type = ct,
                gene      = row$gene,
                TF        = row$TF,
                peak      = row$peak,
                access    = access,
                expr      = expr,
                tf_expr   = tf_expr
            )
        }
        rm(mats); gc()
    }
    bind_rows(plot_data_list) |>
        mutate(
            trio_id     = paste(cell_type, gene, TF, peak, sep = "|"),
            tf_expr_log = log2(tf_expr + 1)
        )
}

make_scatter_pdf = function(plot_df, pdf_path) {
    plot_list = list()
    for (ct in unique(plot_df$cell_type)) {
        ct_df = plot_df |> filter(cell_type == ct)

        ct_plots = list()
        for (tid in unique(ct_df$trio_id)) {
            #   Peak-gene relationship in TRIPOD is only assessed at sufficient
            #   TF expression levels
            trio_data = ct_df |> filter(trio_id == tid, tf_expr > 0)

            ct_plots[[length(ct_plots) + 1]] = ggplot(
                    trio_data,
                    aes(x = access, y = expr, color = tf_expr_log)
                ) +
                geom_point(size = 0.15) +
                scale_color_viridis_c() +
                labs(
                    title = sprintf(
                        "%s\n%s | %s", ct,
                        trio_data$gene[1], trio_data$TF[1]
                    ),
                    x = "Access.", y = "Expr."
                ) +
                theme_bw(base_size = 6) +
                theme(
                    plot.title = element_text(size = 2.5, lineheight = 1.1),
                    legend.position = "none"
                )
        }

        # Pad with blank panels so each cell type occupies a full row
        while (length(ct_plots) < num_examples) {
            ct_plots[[length(ct_plots) + 1]] = ggplot() + theme_void()
        }

        plot_list = c(plot_list, ct_plots)
    }

    n_cell_types = length(unique(plot_df$cell_type))
    pdf(pdf_path, width = 5, height = n_cell_types)
    print(plot_grid(plotlist = plot_list, ncol = num_examples))
    dev.off()
}

################################################################################
#   Main
################################################################################

#-------------------------------------------------------------------------------
#   1. Read trios at FDR < 0.05 (main) and FDR < 0.1 (uniqueness background)
#-------------------------------------------------------------------------------

trio_df_list = list()
background_df_list = list()
for (cell_type in cell_types) {
    this_path = sprintf(trio_paths, cell_type)
    if (file.exists(this_path)) {
        raw = read_parquet_duckdb(this_path, prudence = 'lavish') |>
            filter(condition_on == 'Yj', adj < intersect_fdr_thres) |>
            select(peak, gene, TF, coef, adj, stringency_level) |>
            mutate(cell_type = cell_type) |>
            collect()
        background_df_list[[cell_type]] = raw
        trio_df_list[[cell_type]] = raw |> filter(adj < fdr_thres)
    } else {
        message(
            sprintf(
                "Trios were not found for cell type %s; skipping", cell_type
            )
        )
    }
}

trio_df = bind_rows(trio_df_list)
background_df = bind_rows(background_df_list)

#-------------------------------------------------------------------------------
#   2. Compute sparsity: for each (cell_type, peak, gene), count metacells
#      where both expression > 0 and accessibility > 0
#-------------------------------------------------------------------------------

# Collect unique peak-gene pairs from both sets; background has pairs not in
# trio_df (FDR 0.05-0.1) that also need sparsity counts
pg_needed = bind_rows(trio_df, background_df) |>
    distinct(cell_type, peak, gene)

sparsity_list = list()
for (ct in unique(pg_needed$cell_type)) {
    mats = qs_read(sprintf(prep_path, ct))$metacell_seur
    ct_pg = pg_needed |> filter(cell_type == ct)

    n_nz = integer(nrow(ct_pg))
    for (i in seq_len(nrow(ct_pg))) {
        expr   = mats$rna[, ct_pg$gene[i], drop = TRUE]
        access = mats$peak[, ct_pg$peak[i], drop = TRUE]
        n_nz[i] = sum(expr > 0 & access > 0)
    }

    sparsity_list[[ct]] = ct_pg |> mutate(n_nonzero = n_nz)
    rm(mats); gc()
}

sparsity_df = bind_rows(sparsity_list)

# Join sparsity counts onto both data frames
trio_df = trio_df |>
    filter(cell_type %in% cell_types_keep) |>
    left_join(sparsity_df, by = c("cell_type", "peak", "gene"))
background_df = background_df |>
    filter(cell_type %in% cell_types_keep) |>
    left_join(sparsity_df, by = c("cell_type", "peak", "gene"))

#-------------------------------------------------------------------------------
#   3. Apply sparsity filter & annotate
#-------------------------------------------------------------------------------

message(
    sprintf(
        "Dropping %d of %d trios (%0.2f%%) with < %d nonzero metacells",
        sum(trio_df$n_nonzero < min_nonzero),
        nrow(trio_df),
        100 * mean(trio_df$n_nonzero < min_nonzero),
        min_nonzero
    )
)
trio_df_sparse = trio_df |> filter(n_nonzero >= min_nonzero)
background_df_sparse = background_df |> filter(n_nonzero >= min_nonzero)

# Build peak-gene-TF background at FDR < intersect_fdr_thres (sparsity-filtered),
# either stringency level, one entry per cell type
background_intersect_sparse = lapply(
    split(background_df_sparse, background_df_sparse$cell_type),
    \(df) unique(paste(df$peak, df$gene, df$TF, sep = "|"))
)

trio_df_sparse = annotate_trios(trio_df_sparse, background_intersect_sparse)

#-------------------------------------------------------------------------------
#   4. Export
#-------------------------------------------------------------------------------

trio_df_sparse |>
    select(-n_nonzero) |>
    compute_parquet(out_path)

#-------------------------------------------------------------------------------
#   5. Scatter plots: random FDR < 0.05 trios before and after sparsity
#      filtering, for each stringency level (4 PDFs)
#-------------------------------------------------------------------------------

scatter_data = list()
for (sl in c(1, 2)) {
    # Before sparsity filtering
    pre_trios = trio_df |>
        filter(stringency_level == sl) |>
        group_by(cell_type) |>
        slice_sample(n = num_examples) |>
        ungroup()
    pre_plot_df = make_scatter_data(pre_trios)
    scatter_data[[sprintf("stringency%d_pre_sparsity", sl)]] = pre_plot_df
    make_scatter_pdf(
        pre_plot_df,
        file.path(
            plot_dir,
            sprintf("scatter_stringency%d_pre_sparsity.pdf", sl)
        )
    )

    # After sparsity filtering
    post_trios = trio_df_sparse |>
        filter(stringency_level == sl) |>
        group_by(cell_type) |>
        slice_sample(n = num_examples) |>
        ungroup()
    post_plot_df = make_scatter_data(post_trios)
    scatter_data[[sprintf("stringency%d_post_sparsity", sl)]] = post_plot_df
    make_scatter_pdf(
        post_plot_df,
        file.path(
            plot_dir,
            sprintf("scatter_stringency%d_post_sparsity.pdf", sl)
        )
    )
}

saveRDS(scatter_data, scatter_data_path)

#-------------------------------------------------------------------------------
#   6. Scatter plots for final intersection trios
#      (is_intersect == TRUE at stringency level 1, up to 5 per cell type)
#-------------------------------------------------------------------------------

intersect_trios = trio_df_sparse |>
    filter(is_intersect, stringency_level == 1) |>
    group_by(cell_type) |>
    slice_sample(n = num_examples) |>
    ungroup()

intersect_plot_df = make_scatter_data(intersect_trios)

make_scatter_pdf(
    intersect_plot_df,
    file.path(plot_dir, "scatter_intersect_trios.pdf")
)

session_info()
  