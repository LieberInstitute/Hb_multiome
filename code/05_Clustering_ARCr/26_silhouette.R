########################################################################
## Compute approximate silhouette widths for selected WNN clustering results
##
## This extracts the silhouette-related logic from 06_jaccard_annotated.R and
## intentionally omits pairwise clustering comparisons and Jaccard logic.
########################################################################

library("Seurat")
library("bluster")
library("dplyr")
library("ggplot2")
library("ggbeeswarm")
library("here")
library("scales")
library("stringr")
library("tidyr")

input_rds_dir <- here("processed-data", "05_Clustering_ARCr", "05_rename_idents", "old")
plot_dir <- here("plots", "05_Clustering_ARCr", "26_silhouette")
plot_base_size <- 16

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

object_names <- c(
  "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r1",
  "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r1",
  "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2",
  "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r2"
)

clustering_labels <- tibble(
  object_name = object_names,
  clustering = str_extract(object_names, "k[34]0_C\\.leiden_lsi_r[12]"),
  clustering_short = recode(
    clustering,
    "k30_C.leiden_lsi_r1" = "k30 r1",
    "k40_C.leiden_lsi_r1" = "k40 r1",
    "k30_C.leiden_lsi_r2" = "k30 r2",
    "k40_C.leiden_lsi_r2" = "k40 r2"
  ),
  rds_path = file.path(input_rds_dir, paste0(object_names, ".rds"))
) |>
  mutate(clustering_short = factor(clustering_short, levels = c("k30 r1", "k40 r1", "k30 r2", "k40 r2")))

if (!all(file.exists(clustering_labels$rds_path))) {
  missing_files <- clustering_labels$rds_path[!file.exists(clustering_labels$rds_path)]
  stop("Missing input RDS file(s):\n", paste(missing_files, collapse = "\n"))
}

clean_cluster_label <- function(x) {
  trimws(gsub("\\(\\d*\\.?\\d*%\\)", "", as.character(x)))
}

is_habenula_cluster <- function(x) {
  str_detect(x, "(^|[._ ])(DD_)?[LM]?Hb|Hb[._ ]?ne")
}

save_plot <- function(plot, filename, width = 10, height = 6) {
  ggsave(file.path(plot_dir, filename), plot, width = width, height = height)
}

plot_subtitle <- function(x, width = 88) {
  stringr::str_wrap(x, width = width)
}

compute_one_silhouette <- function(rds_path, clustering, clustering_short) {
  message("Processing ", clustering_short)
  seurat_obj <- readRDS(rds_path)

  clusters <- clean_cluster_label(Idents(seurat_obj))
  embeddings <- Embeddings(seurat_obj, reduction = "integrated.harmony")
  donors <- seurat_obj[["orig.ident"]][, 1]
  cells <- colnames(seurat_obj)
  hb_keep <- is_habenula_cluster(clusters)

  all_sil <- approxSilhouette(embeddings, clusters) |>
    as.data.frame() |>
    mutate(
      clustering = clustering,
      clustering_short = clustering_short,
      scope = "all_cells",
      cell = cells,
      donor = donors,
      cluster = clusters,
      is_habenula = hb_keep,
      other = clean_cluster_label(other),
      closest = if_else(width > 0, cluster, other)
    ) |>
    select(clustering, clustering_short, scope, cell, donor, cluster, is_habenula, other, closest, width)

  hb_sil <- approxSilhouette(embeddings[hb_keep, , drop = FALSE], clusters[hb_keep]) |>
    as.data.frame() |>
    mutate(
      clustering = clustering,
      clustering_short = clustering_short,
      scope = "habenula_only",
      cell = cells[hb_keep],
      donor = donors[hb_keep],
      cluster = clusters[hb_keep],
      is_habenula = TRUE,
      other = clean_cluster_label(other),
      closest = if_else(width > 0, cluster, other)
    ) |>
    select(clustering, clustering_short, scope, cell, donor, cluster, is_habenula, other, closest, width)

  donor_counts <- tibble(
    clustering = clustering,
    clustering_short = clustering_short,
    cluster = clusters,
    donor = donors
  ) |>
    count(clustering, clustering_short, cluster, donor, name = "n")

  summary <- tibble(
    clustering = clustering,
    clustering_short = clustering_short,
    n_cells = length(cells),
    n_clusters = n_distinct(clusters),
    n_hb_clusters = n_distinct(clusters[hb_keep]),
    n_donors = n_distinct(donors)
  )

  list(
    silhouette_data = bind_rows(all_sil, hb_sil),
    donor_counts = donor_counts,
    summary = summary
  )
}

plot_silhouette_swarm <- function(data, scope_value, filename) {
  plt <- data |>
    filter(clustering == "k30_C.leiden_lsi_r2", scope == scope_value) |>
    mutate(cluster = factor(cluster, levels = sort(unique(cluster)))) |>
    ggplot(aes(x = cluster, y = width, colour = closest)) +
    ggbeeswarm::geom_quasirandom(method = "smiley", alpha = 0.4) +
    labs(
      title = paste("Approximate silhouette widths:", "k30 r2", scope_value),
      subtitle = plot_subtitle("Original diagnostic from earlier scripts. Useful QC, but not the most convincing comparison plot because it only shows k30 r2."),
      x = NULL,
      y = "Approximate silhouette width"
    ) +
    theme_minimal(base_size = plot_base_size) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "none"
    )

  save_plot(plt, filename, height = 7, width = 12)
  plt
}

silhouette_results <- vector("list", nrow(clustering_labels))

for (i in seq_len(nrow(clustering_labels))) {
  silhouette_results[[i]] <- compute_one_silhouette(
    clustering_labels$rds_path[i],
    clustering_labels$clustering[i],
    clustering_labels$clustering_short[i]
  )
  gc()
}

silhouette_data <- bind_rows(lapply(silhouette_results, `[[`, "silhouette_data")) |>
  mutate(clustering_short = factor(clustering_short, levels = c("k30 r1", "k40 r1", "k30 r2", "k40 r2")))
donor_counts <- bind_rows(lapply(silhouette_results, `[[`, "donor_counts")) |>
  mutate(clustering_short = factor(clustering_short, levels = c("k30 r1", "k40 r1", "k30 r2", "k40 r2")))
run_summary <- bind_rows(lapply(silhouette_results, `[[`, "summary")) |>
  mutate(clustering_short = factor(clustering_short, levels = c("k30 r1", "k40 r1", "k30 r2", "k40 r2")))

rm(silhouette_results)
gc()

########################################################################
## Stable outputs
########################################################################

plot_silhouette_swarm(
  silhouette_data,
  "all_cells",
  "k30_C.leiden_lsi_r2_silhouette_all_cells.png"
)
plot_silhouette_swarm(
  silhouette_data,
  "habenula_only",
  "k30_C.leiden_lsi_r2_silhouette_habenula_only.png"
)

hb_cluster_totals <- silhouette_data |>
  filter(scope == "all_cells", is_habenula) |>
  group_by(clustering, clustering_short) |>
  summarise(total_hb_clusters = n_distinct(cluster), .groups = "drop")

threshold_summary <- silhouette_data |>
  filter(scope == "all_cells", is_habenula) |>
  tidyr::crossing(threshold = seq(-0.25, 0.75, by = 0.025)) |>
  group_by(clustering, clustering_short, threshold) |>
  summarise(n_hb_clusters = n_distinct(cluster[width >= threshold]), .groups = "drop") |>
  left_join(hb_cluster_totals, by = c("clustering", "clustering_short")) |>
  mutate(prop_hb_clusters = n_hb_clusters / total_hb_clusters)

plt_hb_threshold <- ggplot(threshold_summary, aes(x = threshold, y = prop_hb_clusters, color = clustering_short)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Hb clusters represented above silhouette threshold",
    subtitle = plot_subtitle("Higher curves retain a larger fraction of Hb clusters as the silhouette threshold increases. Convincing only if k30 r2 separates from alternatives over a broad threshold range."),
    x = "Approximate silhouette width threshold",
    y = "Proportion of Hb clusters",
    color = "Clustering"
  ) +
  theme_minimal(base_size = plot_base_size) +
  theme(legend.position = "bottom")
save_plot(plt_hb_threshold, "hb_cluster_proportion_by_silhouette_threshold.pdf")

donor_purity <- donor_counts |>
  group_by(clustering, clustering_short, cluster) |>
  mutate(
    cluster_n = sum(n),
    donor_prop = n / cluster_n,
    n_donors = n_distinct(donor),
    max_donor_prop = max(donor_prop),
    inverse_simpson = 1 / sum(donor_prop^2),
    normalized_entropy = -sum(donor_prop * log(donor_prop)) / log(n_distinct(donor_counts$donor))
  ) |>
  ungroup() |>
  select(
    clustering,
    clustering_short,
    cluster,
    cluster_n,
    n_donors,
    max_donor_prop,
    inverse_simpson,
    normalized_entropy
  ) |>
  distinct()

plt_donor_balance <- ggplot(donor_purity, aes(
  x = max_donor_prop,
  y = normalized_entropy,
  size = cluster_n,
  color = clustering_short
)) +
  geom_point(alpha = 0.7) +
  labs(
    title = "Cluster donor balance",
    subtitle = plot_subtitle("Lower x and higher y indicate less donor-dominated clusters. Useful diagnostic, but not especially decisive for favoring k30 r2."),
    x = "Largest donor proportion in cluster",
    y = "Normalized donor entropy",
    size = "Cells in cluster",
    color = "Clustering"
  ) +
  theme_minimal(base_size = plot_base_size) +
  theme(legend.position = "bottom")
save_plot(plt_donor_balance, "cluster_donor_balance.pdf")

plt_cluster_max_donor <- ggplot(donor_purity, aes(x = clustering_short, y = max_donor_prop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(size = cluster_n), width = 0.15, alpha = 0.65) +
  labs(
    title = "Maximum donor contribution per cluster",
    subtitle = plot_subtitle("Lower values indicate less donor-dominated clusters. Simple and interpretable, but cluster weighting can overemphasize small clusters."),
    x = NULL,
    y = "Largest donor proportion in cluster",
    size = "Cells in cluster"
  ) +
  theme_minimal(base_size = plot_base_size) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(plt_cluster_max_donor, "cluster_max_donor_proportion_boxplot.pdf")

########################################################################
## Experimental and exploratory plots
########################################################################

silhouette_cluster_summary <- silhouette_data |>
  mutate(scope_label = case_when(
    scope == "all_cells" & is_habenula ~ "Hb clusters vs all cells",
    scope == "habenula_only" ~ "Hb clusters vs Hb only",
    scope == "all_cells" ~ "All clusters vs all cells"
  )) |>
  group_by(clustering, clustering_short, scope_label, cluster) |>
  summarise(
    n_cells = n(),
    median_width = median(width),
    mean_width = mean(width),
    prop_positive = mean(width > 0),
    prop_gt_0.05 = mean(width > 0.05),
    prop_gt_0.10 = mean(width > 0.10),
    .groups = "drop"
  )

candidate_cluster_table <- silhouette_cluster_summary |>
  filter(scope_label == "Hb clusters vs all cells") |>
  select(clustering, clustering_short, cluster, n_cells, median_width, prop_positive, prop_gt_0.10) |>
  left_join(
    donor_purity |>
      transmute(clustering, cluster, max_donor_prop, normalized_entropy),
    by = c("clustering", "cluster")
  ) |>
  mutate(joint_good = median_width > 0.10 & max_donor_prop <= 0.35)

plt_joint_scatter <- ggplot(candidate_cluster_table, aes(x = median_width, y = max_donor_prop)) +
  geom_rect(
    xmin = 0.10,
    xmax = Inf,
    ymin = -Inf,
    ymax = 0.35,
    fill = "grey85",
    alpha = 0.5
  ) +
  geom_vline(xintercept = 0.10, linetype = "dashed") +
  geom_hline(yintercept = 0.35, linetype = "dashed") +
  geom_point(aes(size = n_cells, color = joint_good), alpha = 0.75) +
  facet_wrap(vars(clustering_short), nrow = 2) +
  labs(
    title = "Experimental: Hb cluster quality screen",
    subtitle = plot_subtitle("Shaded area: median silhouette > 0.10 and largest donor proportion <= 0.35. Initially promising for k30 r2, but the donor cutoff is arbitrary and should not be the lead evidence."),
    x = "Median silhouette width",
    y = "Largest donor proportion",
    size = "Cells",
    color = "Pass"
  ) +
  theme_minimal(base_size = plot_base_size) +
  theme(legend.position = "bottom")
save_plot(plt_joint_scatter, "experimental_hb_joint_quality_scatter.pdf", width = 11, height = 8)

joint_pass_summary <- candidate_cluster_table |>
  group_by(clustering_short) |>
  summarise(
    n_hb_clusters = n(),
    n_pass = sum(joint_good),
    prop_pass = mean(joint_good),
    .groups = "drop"
  )

plt_joint_pass <- ggplot(joint_pass_summary, aes(x = clustering_short, y = prop_pass)) +
  geom_col() +
  geom_text(aes(label = paste0(n_pass, "/", n_hb_clusters)), vjust = -0.4, size = 5) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 0.5)) +
  labs(
    title = "Experimental: Hb clusters passing joint quality screen",
    subtitle = plot_subtitle("k30 r2 has the largest fraction passing median silhouette > 0.10 and max donor proportion <= 0.35. Visually decisive, but too dependent on arbitrary cutoffs to rely on alone."),
    x = NULL,
    y = "Hb clusters passing"
  ) +
  theme_minimal(base_size = plot_base_size)
save_plot(plt_joint_pass, "experimental_hb_joint_quality_pass_barplot.pdf")

sensitivity_grid <- tidyr::crossing(
  silhouette_cutoff = seq(0.00, 0.20, by = 0.025),
  donor_cutoff = seq(0.25, 0.50, by = 0.025)
)

sensitivity_summary <- candidate_cluster_table |>
  tidyr::crossing(sensitivity_grid) |>
  mutate(pass = median_width > silhouette_cutoff & max_donor_prop <= donor_cutoff) |>
  group_by(clustering_short, silhouette_cutoff, donor_cutoff) |>
  summarise(prop_pass = mean(pass), .groups = "drop")

plt_sensitivity <- sensitivity_summary |>
  filter(donor_cutoff %in% c(0.30, 0.35, 0.40, 0.45)) |>
  ggplot(aes(x = silhouette_cutoff, y = prop_pass, color = clustering_short)) +
  geom_line(linewidth = 1) +
  facet_wrap(vars(donor_cutoff), labeller = label_both) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Experimental: sensitivity of Hb joint quality screen",
    subtitle = plot_subtitle("k30 r2 looks best near stricter donor cutoffs, but the advantage weakens as donor thresholds relax. Useful as a transparency check, not the cleanest final figure."),
    x = "Median silhouette cutoff",
    y = "Hb clusters passing",
    color = "Clustering"
  ) +
  theme_minimal(base_size = plot_base_size) +
  theme(legend.position = "bottom")
save_plot(plt_sensitivity, "experimental_hb_joint_quality_sensitivity.pdf", width = 11, height = 8)

cluster_donor_metrics <- donor_counts |>
  group_by(clustering, clustering_short, cluster) |>
  mutate(
    cluster_n = sum(n),
    donor_prop = n / cluster_n,
    max_donor_prop = max(donor_prop),
    inverse_simpson = 1 / sum(donor_prop^2),
    normalized_entropy = -sum(donor_prop * log(donor_prop)) / log(n_distinct(donor_counts$donor))
  ) |>
  ungroup()

hb_cell_scores <- silhouette_data |>
  filter(scope == "all_cells", is_habenula) |>
  left_join(
    cluster_donor_metrics |>
      select(clustering, cluster, donor, donor_prop, max_donor_prop, inverse_simpson, normalized_entropy),
    by = c("clustering", "cluster", "donor")
  ) |>
  mutate(clustering_short = factor(clustering_short, levels = c("k30 r1", "k40 r1", "k30 r2", "k40 r2")))

plt_hb_cell_ecdf <- ggplot(hb_cell_scores, aes(x = max_donor_prop, color = clustering_short)) +
  stat_ecdf(linewidth = 1.2) +
  labs(
    title = "Experimental: Hb-cell exposure to donor-dominated clusters",
    subtitle = plot_subtitle("Among most proportions p, k30 r2 has a higher proportion of habenula cells present in clusters satisfying max donor proportion < p. Likely one of the more convincing visualizations."),
    x = "Largest donor proportion in assigned cluster",
    y = "Fraction of Hb cells",
    color = "Clustering"
  ) +
  theme_minimal(base_size = plot_base_size) +
  theme(legend.position = "bottom")
save_plot(plt_hb_cell_ecdf, "experimental_hb_cell_donor_dominance_ecdf.pdf")

plt_hb_cell_violin <- ggplot(hb_cell_scores, aes(x = clustering_short, y = max_donor_prop)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.35) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.75) +
  labs(
    title = "Experimental: cell-weighted donor dominance of Hb clusters",
    subtitle = plot_subtitle("Each Hb cell contributes its assigned cluster's max donor proportion. k30 r2 has a favorable central distribution, though this is less visually decisive than the ECDF."),
    x = NULL,
    y = "Largest donor proportion in assigned cluster"
  ) +
  theme_minimal(base_size = plot_base_size)
save_plot(plt_hb_cell_violin, "experimental_hb_cell_donor_dominance_violin.pdf")

message("Process completed. Plots written to:\n", plot_dir)
run_summary
