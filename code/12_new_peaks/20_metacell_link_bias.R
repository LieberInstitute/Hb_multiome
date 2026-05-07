library(tidyverse)
library(Seurat)
library(Signac)
library(here)
library(qs2)
library(sessioninfo)
library(duckplyr)
library(matrixStats)

seur_path = here(
    "processed-data", "11_link_prep", "05_metacell_aggregate",
    "meta_seur.qs2"
)
link_path = here(
    'processed-data', '12_new_peaks', '17_filter_links',
    'metacell_filtered_data.parquet'
)
out_path = here(
    'processed-data', '12_new_peaks', '20_metacell_link_bias',
    'metacell_robustness_scores.parquet'
)
plot_dir = here("plots", "12_new_peaks", "20_metacell_link_bias")

atac_assay = "ATAC"
rna_assay  = "RNA"
#   Fraction of metacells to sample per permutation. At 50%, even MHb.3 (10
#   metacells) gets 5 per draw — enough for a meaningful correlation — while
#   still testing robustness to dropping half the data.
subsample_frac = 0.5
#   Number of subsampling permutations
n_perm = 50

num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(num_cores)) num_cores = 1L
duckplyr::db_exec(sprintf("SET threads = %d", num_cores))
fallback_config(info = FALSE)

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
set.seed(0)

seur = qs_read(seur_path)

link_df = read_parquet_duckdb(link_path) |>
    collect()

stopifnot(all(link_df$peak %in% rownames(seur[[atac_assay]])))
stopifnot(all(link_df$gene %in% rownames(seur[[rna_assay]])))

#   Z-score rows of a matrix across columns (metacells); rows with near-zero sd
#   are set to NA so they contribute NA correlations rather than NaN/Inf
zscore_rows = function(mat) {
    mu  = rowMeans(mat)
    sds = rowSds(mat)
    sds[sds < 1e-10] = NA_real_
    (mat - mu) / sds
}

ct_meta = seur@meta.data |>
    as_tibble(rownames = "metacell")

score_df_list = vector("list", length(unique(link_df$cell_type)))
names(score_df_list) = unique(link_df$cell_type)

for (this_cell_type in unique(link_df$cell_type)) {
    message(sprintf("Processing %s ...", this_cell_type))

    ct_cells = ct_meta |>
        filter(refined_mid_cluster == this_cell_type) |>
        pull(metacell)
    ct_links = link_df |>
        filter(cell_type == this_cell_type)

    n_cells   = length(ct_cells)
    samp_size = floor(subsample_frac * n_cells)

    #   Pull expression matrices (links x metacells) and z-score across metacells.
    #   z-scored correlations reduce the confound of expression magnitude, which
    #   is relevant since metacell size varies.
    peak_mat = as.matrix(
        GetAssayData(seur, assay = atac_assay, layer = "data")[ct_links$peak, ct_cells]
    )
    gene_mat = as.matrix(
        GetAssayData(seur, assay = rna_assay, layer = "data")[ct_links$gene, ct_cells]
    )
    peak_z = zscore_rows(peak_mat)
    gene_z = zscore_rows(gene_mat)

    #   Full-sample correlation for each peak-gene pair (Pearson via dot product
    #   of z-scores, divided by n)
    full_cor = rowMeans(peak_z * gene_z, na.rm = TRUE)

    #   Permutation loop: each iteration subsamples half the metacells and
    #   recomputes correlations, building a distribution to characterize
    #   how stable the full-sample correlation is.
    perm_cors = matrix(NA_real_, nrow = nrow(peak_z), ncol = n_perm)
    for (p in seq_len(n_perm)) {
        idx = sample(n_cells, samp_size)
        perm_cors[, p] = rowMeans(peak_z[, idx] * gene_z[, idx], na.rm = TRUE)
    }

    #   Robustness metrics — both avoid the ratio instability near zero that
    #   plagues the donor_score approach:
    #
    #   sign_flip_rate: fraction of permutations where the subsampled correlation
    #     has the opposite sign to the full correlation. Bounded [0, 1]. Values
    #     near 0.5 indicate no reliable directional signal.
    #
    #   instability: mean |full_cor - perm_cor| across permutations. In the same
    #     units as a correlation coefficient, so directly interpretable.
    fc_mat         = matrix(full_cor, nrow = length(full_cor), ncol = n_perm)
    sign_flip_rate = rowMeans(fc_mat * perm_cors <= 0, na.rm = TRUE)
    instability    = rowMeans(abs(fc_mat - perm_cors), na.rm = TRUE)

    score_df_list[[this_cell_type]] = tibble(
        cell_type      = this_cell_type,
        peak           = ct_links$peak,
        gene           = ct_links$gene,
        full_cor       = full_cor,
        sign_flip_rate = sign_flip_rate,
        instability    = instability
    )
}

score_df = bind_rows(score_df_list)

#   Join back original link metadata
score_df = score_df |>
    left_join(link_df, by = c("cell_type", "peak", "gene"))

score_df |>
    compute_parquet(out_path)

session_info()
