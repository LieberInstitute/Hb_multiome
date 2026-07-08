#Quantifying the RNAscope data of the VGAT probes in the habenula


library(ggplot2)
library(here)


plot_path = here('plots','14_vgat_rnascope', '01_halo_data')
if (!dir.exists(plot_path)) dir.create(plot_path)
new_data_path = here('processed-data','14_vgat_rnascope', '01_halo_data')
if (!dir.exists(new_data_path)) dir.create(new_data_path)

halo_path = here('processed-data','14_vgat_rnascope','HALO_Trimmed')


all_files = list.files(halo_path)
csv_files = all_files[grepl('.csv',all_files)]


# Loop through all samples, collecting threshold sweep data
threshold_list <- list()

for (sample_id in seq_along(csv_files)) {

  current_data = data.table::fread(file.path(halo_path, csv_files[sample_id]))
  sample_name = strsplit(csv_files[sample_id], split = '_', fixed = TRUE)[[1]][1]
  cat(sprintf("\n=== Processing sample %d: %s ===\n", sample_id, sample_name))

  # Compute centroids
  current_data[, x_centroid := (XMin + XMax) / 2]
  current_data[, y_centroid := (YMin + YMax) / 2]


  # Classify habenula (POU4F1 >= 10) and VGAT+ (SLC32A1 >= 1)
  current_data[, pou4f1_pos := `POU4F1_690 (Opal 690) Copies` >= 10]
  current_data[, vgat_pos := `SLC32A1_570 (Opal 570) Copies` >= 1]

  # Combined annotation
  current_data[, annotation := factor(
    data.table::fcase(
      pou4f1_pos & vgat_pos, "POU4F1+ and VGAT+",
      pou4f1_pos & !vgat_pos, "POU4F1+ only",
      !pou4f1_pos & vgat_pos, "VGAT+ only",
      !pou4f1_pos & !vgat_pos, "Neither"
    ),
    levels = c("POU4F1+ and VGAT+", "POU4F1+ only", "VGAT+ only", "Neither")
  )]

  annotation_colors <- c(
    "POU4F1+ and VGAT+" = "#020202",
    "POU4F1+ only" = "#a80fe9",
    "VGAT+ only" = "#f12828",
    "Neither" = "#e0dede"
  )

  # Plot all objects colored by combined annotation
  p1 = ggplot(current_data, aes(x = x_centroid, y = y_centroid, color = annotation)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = annotation_colors) +
    scale_y_reverse() +
    labs(title = sprintf("POU4F1 & VGAT classification - %s", sample_name),
        x = "X position", y = "Y position", color = NULL) +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    theme(
      aspect.ratio = 1,
      legend.text = element_text(size = 12),
      legend.key.size = unit(1, "cm")
    )



  # Define habenula bounding box using quantiles to exclude outliers
  hb_objects <- current_data[pou4f1_pos == TRUE]

  # Use 5th and 95th percentiles to trim outliers
  hb_x_range <- quantile(hb_objects$x_centroid, probs = c(0.05, 0.95))
  hb_y_range <- quantile(hb_objects$y_centroid, probs = c(0.05, 0.95))

  cat("Habenula bounding box (5%-95% quantiles):\n")
  cat("  X range:", hb_x_range[1], "-", hb_x_range[2], "\n")
  cat("  Y range:", hb_y_range[1], "-", hb_y_range[2], "\n")


  # Moving window analysis
  # Window spans the full x-range of the habenula, moves downward in y
  # Start at the top of the habenula (min y) and scan down to the bottom of the tissue

  window_x_min <- hb_x_range[1] 
  window_x_max <- hb_x_range[2] 

  window_height <- 100
  step_size <- 25

  # Y range: start at top of habenula, go to bottom of tissue
  y_start <- hb_y_range[1]
  y_end <- max(current_data$y_centroid)

  # Generate window positions
  window_starts <- seq(y_start, y_end - window_height, by = step_size)



  # Bounding box plot with habenula boundary and first window annotation
  first_window_ymin <- window_starts[1]
  first_window_ymax <- window_starts[1] + window_height

  p2 = ggplot(current_data, aes(x = x_centroid, y = y_centroid, color = pou4f1_pos)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = c("FALSE" = "#e0dede", "TRUE" = "#a80fe9"),
                      labels = c("Non-POU4F1", "POU4F1 (copies >= 10)")) +
    annotate("rect",
            xmin = hb_x_range[1], xmax = hb_x_range[2],
            ymin = hb_y_range[1], ymax = hb_y_range[2],
            fill = NA, color = "black", linewidth = 0.8, linetype = "dashed") +
    annotate("rect",
            xmin = window_x_min, xmax = window_x_max,
            ymin = first_window_ymin, ymax = first_window_ymax,
            fill = "#4472C4", alpha = 0.3, color = "#4472C4", linewidth = 0.6) +
    annotate("text", x = window_x_max + 100, y = (first_window_ymin + first_window_ymax) / 2,
            label = "First window", hjust = 0, size = 5, color = "#4472C4") +
    geom_hline(yintercept = hb_y_range[2], linetype = "dashed", color = "grey40", size = 1) +
    annotate("text", x = max(current_data$x_centroid) * 0.85, y = hb_y_range[2],
            label = "Hb boundary", vjust = -0.5, size = 5, color = "grey40") +
    scale_y_reverse() +
    labs(title = paste0("Habenula bounding box (5%-95% quantiles) - ", sample_name),
        x = "X position", y = "Y position", color = NULL) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.key.size = unit(1, "cm")
    )




  # Recompute moving window with mean and sd for each probe
  results <- data.table::rbindlist(lapply(window_starts, function(y_top) {
    in_window <- current_data[
      x_centroid >= window_x_min & x_centroid <= window_x_max &
      y_centroid >= y_top & y_centroid < y_top + window_height
    ]
    
    data.table::data.table(
      y_center = y_top + window_height / 2,
      n_objects = nrow(in_window),
      SLC32A1_mean = mean(in_window$`SLC32A1_570 (Opal 570) Copies`, na.rm = TRUE),
      SLC32A1_sd = sd(in_window$`SLC32A1_570 (Opal 570) Copies`, na.rm = TRUE),
      POU4F1_mean = mean(in_window$`POU4F1_690 (Opal 690) Copies`, na.rm = TRUE),
      POU4F1_sd = sd(in_window$`POU4F1_690 (Opal 690) Copies`, na.rm = TRUE),
      MBP_mean = mean(in_window$`MBP_520 (Opal 520) Copies`, na.rm = TRUE),
      MBP_sd = sd(in_window$`MBP_520 (Opal 520) Copies`, na.rm = TRUE)
    )
  }))

  # Pivot to long format
  results_long <- data.table::melt(results, 
    id.vars = c("y_center", "n_objects"),
    measure.vars = patterns("_mean$", "_sd$"),
    variable.name = "probe",
    value.name = c("mean_copies", "sd_copies")
  )
  results_long[, probe := c("SLC32A1", "POU4F1", "MBP")[probe]]


  p3 = ggplot(results_long, aes(x = y_center, y = mean_copies, color = probe, fill = probe)) +
    geom_ribbon(aes(ymin = pmax(0, mean_copies - sd_copies), ymax = mean_copies + sd_copies), 
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = hb_y_range[2], linetype = "dashed", color = "grey40", size = 1) +
    annotate("text", x = hb_y_range[2], y = max(results_long$mean_copies, na.rm = TRUE) * 0.95,
            label = "Hb boundary", hjust = -0.1, size = 5, color = "grey40") +
    scale_color_manual(values = c("SLC32A1" = "#f12828", "SLC17A6" = "#D95F02", "POU4F1" = "#a80fe9", "MBP" = "#1B9E77")) +
    scale_fill_manual(values = c("SLC32A1" = "#f12828", "SLC17A6" = "#D95F02", "POU4F1" = "#a80fe9", "MBP" = "#1B9E77")) +
    labs(title = sprintf("Moving window: mean probe copies ± SD - %s", sample_name),
        x = "Y position (top → bottom)",
        y = "Mean probe copies per cell",
        color = "Probe", fill = "Probe") +
    theme_minimal()



  # Classify objects as inside/outside the habenula bounding box
  current_data[, in_hb_region := x_centroid >= hb_x_range[1] & x_centroid <= hb_x_range[2] &
                y_centroid >= hb_y_range[1] & y_centroid <= hb_y_range[2]]


  # Barplot: proportion of VGAT+ and POU4F1+ objects inside vs outside Hb region
  count_data <- current_data[, .(
    `VGAT+` = sum(vgat_pos) / .N,
    `POU4F1+` = sum(pou4f1_pos) / .N
  ), by = .(region = ifelse(in_hb_region, "Inside Hb", "Outside Hb"))]

  count_long <- data.table::melt(count_data, id.vars = "region",
                                variable.name = "category", value.name = "proportion")

  p4 = ggplot(count_long, aes(x = category, y = proportion, fill = region)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Inside Hb" = "#D73027", "Outside Hb" = "grey70")) +
    labs(title = sprintf("VGAT+ and POU4F1+ proportion by region - %s", sample_name),
        x = NULL, y = "Proportion of cells", fill = NULL) +
    theme_minimal() +
    theme(legend.text = element_text(size = 12),
          legend.key.size = unit(0.8, "cm"))


  # Proportion of Hb region objects that are VGAT+ across varying thresholds
  hb_data <- current_data[in_hb_region == TRUE]
  thresholds <- seq(1, 50, by = 1)

  threshold_results <- data.table::data.table(
    threshold = thresholds,
    proportion_vgat = sapply(thresholds, function(t) {
      sum(hb_data$`SLC32A1_570 (Opal 570) Copies` >= t) / nrow(hb_data)
    }),
    sample = sample_name
  )

  p5 = ggplot(threshold_results, aes(x = threshold, y = proportion_vgat)) +
    geom_line(linewidth = 1, color = "#f12828") +
    geom_point(size = 0.8, color = "#f12828") +
    labs(title = sprintf("Proportion of Hb cells that are VGAT+ by threshold - %s", sample_name),
         x = "SLC32A1 copy threshold (>=)",
         y = "Proportion of habenula cells") +
    theme_minimal()


  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)


  ggsave(plot = p1, path = plot_path, filename = sprintf("%s_combined_annotation.png", sample_name), 
  width = 8, height = 8, device = 'png')

  ggsave(plot = p2, path = plot_path, filename = sprintf("%s_habenula_box.png", sample_name), 
  width = 8, height = 8, device = 'png')

  ggsave(plot = p3, path = plot_path, filename = sprintf("%s_moving_window_averages.pdf", sample_name), 
  width = 6, height = 4, device = 'pdf')

  ggsave(plot = p4, path = plot_path, filename = sprintf("%s_object_count_in_habenula_barplot.pdf", sample_name), 
  width = 6, height = 4, device = 'pdf')

  ggsave(plot = p5, path = plot_path, filename = sprintf("%s_vgat_threshold_sweep.pdf", sample_name), 
  width = 6, height = 4, device = 'pdf')

  # Collect threshold data for summary plot
  threshold_list[[sample_id]] <- threshold_results
}


# Summary plot across all donors
all_thresholds <- data.table::rbindlist(threshold_list)

p6 = ggplot(all_thresholds, aes(x = threshold, y = proportion_vgat, color = sample)) +
  geom_line(linewidth = 1) +
  geom_point(size = 0.8) +
  labs(title = "Proportion of Hb cells that are VGAT+ by threshold - All donors",
       x = "SLC32A1 copy threshold (>=)",
       y = "Proportion of habenula cells",
       color = "Donor") +
  theme_minimal() +
  theme(legend.text = element_text(size = 11))

print(p6)

ggsave(plot = p6, path = plot_path, filename = "all_donors_vgat_threshold_sweep.pdf",
width = 7, height = 5, device = 'pdf')
