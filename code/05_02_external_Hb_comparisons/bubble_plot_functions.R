
#Found I was reusing the same plotting code over multiple scripts, just consolidate here




#Custom bubble plot, gets mean expression per cluster for a gene, plots the z-score of that across the clusters
get_bubble_plot = function(seurat_object, top_markers, sample_name, group_col = "meta_cluster"){
  # Extract expression data and metadata
  expr_data <- FetchData(seurat_object, vars = top_markers, layer = "data")
  metadata <- seurat_object@meta.data

  # Combine into a data frame
  plot_data <- cbind(expr_data, group_var = metadata[[group_col]]) %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = -group_var, names_to = "gene", values_to = "expression")

  # Calculate mean expression and percent expressing per cluster
  summary_data <- plot_data %>% filter(group_var != 'outliers') %>%
    group_by(gene, group_var) %>%
    summarise(
      mean_expression = mean(expression),
      pct_expressing = sum(expression > 0) / n() * 100,
      .groups = "drop"
    ) %>%
    # Calculate z-score of mean_expression per gene across clusters
    group_by(gene) %>%
    mutate(mean_expression_zscore = scale(mean_expression)[,1]) %>%
    ungroup()

  # Set factor levels to control axis order
  summary_data$gene <- factor(summary_data$gene, levels = top_markers)
  summary_data$group_var <- factor(summary_data$group_var, 
                                      levels = sort(unique(summary_data$group_var)))

  # Create bubble plot
  p1 = ggplot(summary_data, aes(x = gene, y = group_var, size = pct_expressing, color = mean_expression)) +
    geom_point() +
    scale_color_gradient2(low = "white", high = "red", name = "Mean Expression") +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() + ggtitle(sample_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = group_col, size = "% Expressing", color = "Mean Expression")
    
  p2 = ggplot(summary_data, aes(x = gene, y = group_var, size = pct_expressing, color = mean_expression_zscore)) +
    geom_point() +
    scale_color_gradient2(low = "blue", mid = 'white', high = "red", name = "Mean Exp. z-score") +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() + ggtitle(sample_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = group_col, size = "% Expressing", color = "Mean Exp. z-score")
  return(list(p1, p2))
}



#Same function but for SCE objects

get_bubble_plot_sce = function(sce_object, top_markers, sample_name, group_col = "meta_cluster", group_order = NULL, exp_assay = 'cpm'){
  # Extract expression data and metadata
  expr_data <- assay(sce_object, exp_assay)[top_markers, ]
  metadata <- colData(sce_object)

  # Convert to data frame for plotting
  plot_data <- as.data.frame(t(as.matrix(expr_data))) %>%
    tibble::rownames_to_column("cell_id") %>%
    cbind(group_var = metadata[[group_col]]) %>%
    tidyr::pivot_longer(cols = -c(cell_id, group_var), names_to = "gene", values_to = "expression")

  # Calculate mean expression and percent expressing per cluster
  summary_data <- plot_data %>% filter(group_var != 'outliers') %>%
    group_by(gene, group_var) %>%
    summarise(
      mean_expression = mean(expression),
      pct_expressing = sum(expression > 0) / n() * 100,
      .groups = "drop"
    ) %>%
    # Calculate z-score of mean_expression per gene across clusters
    group_by(gene) %>%
    mutate(mean_expression_zscore = scale(mean_expression)[,1]) %>%
    ungroup()

  # Set factor levels to control axis order
  summary_data$gene <- factor(summary_data$gene, levels = top_markers)
  if(is.null(group_order)){
    summary_data$group_var <- factor(summary_data$group_var, 
                                      levels = sort(unique(summary_data$group_var)))
  } else {
    summary_data$group_var <- factor(summary_data$group_var, levels = group_order)
  }
  max_val = max(summary_data$mean_expression)
  # Create bubble plot
  p1 = ggplot(summary_data, aes(x = gene, y = group_var, size = pct_expressing, color = mean_expression)) +
    geom_point() +
    scale_color_gradientn(colours = c("grey70", viridisLite::viridis(256)[20:256]),
      limits = c(0, max_val)) +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() + ggtitle(sample_name) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Gene", y = group_col, size = "% Expression", color = "Mean Expression")

  p2 = ggplot(summary_data, aes(x = gene, y = group_var, size = pct_expressing, color = mean_expression_zscore)) +
    geom_point() +
    scale_color_gradient2(low = "blue", mid = 'white', high = "red", midpoint = 0,
    #limits = c(-2, 2),
    #oob = scales::squish,
    name = "Mean Exp. z-score") +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() + ggtitle(sample_name) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Gene", y = group_col, size = "% Expressing", color = "Mean Exp. z-score")
  return(list(p1, p2))
}




