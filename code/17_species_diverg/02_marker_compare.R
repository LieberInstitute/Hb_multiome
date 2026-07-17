#Comparing marker stats between human and mouse for the consensus cell-type annotations

library(ggplot2)
library(ggrepel)
library(dplyr)
library(here)
library(MetaMarkers)

here::here()


#Path to save any generated data
new_data_path = here('processed-data', '17_species_diverg', '02_marker_compare')
#Path to plot directory
plot_path = here('plots','17_species_diverg', '02_marker_compare')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)

#colors
source(here('code','05_03_annotation_adjustments','celltype_colors.R'))

human_marker_path = here('processed-data', '05_03_annotation_adjustments', '08_metamarkers')
mouse_marker_path = here('processed-data','17_species_diverg','01_hashikawa_markers')


#human_marker_stats = readRDS(file = paste0(human_marker_path, '/marker_stats_combo.rds'))
#mouse_marker_stats = readRDS(file = paste0(mouse_marker_path, '/marker_stats_combo.rds'))

human_markers = read_meta_markers(paste0(human_marker_path, '/multiome_no_thal_mid_meta_markers.csv.gz'))
mouse_markers = read_meta_markers(paste0(mouse_marker_path, '/mouse_meta_markers.csv.gz'))

#Filter for shared genes
shared_genes = intersect(unique(human_markers$gene), unique(mouse_markers$gene))
length(shared_genes)

human_markers = human_markers[human_markers$gene %in% shared_genes, ]
mouse_markers = mouse_markers[mouse_markers$gene %in% shared_genes, ]


human_markers |> filter(rank <= 75) |> View()
mouse_markers |> filter(rank <= 75) |> View()


#Update the mouse cell-type names
mouse_markers$cell_type = recode(mouse_markers$cell_type, 
  'MHb.1' = 'MHb_A', 'MHb.2' = 'MHb_B', 'LHb.2.7' = 'LHb_A', 
  'LHb.1.3.4' = 'LHb_B', 'LHb.4' = 'LHb_C')


plot_de_stats = function(hu_markers, mou_markers, stat_of_interest, 
  human_celltype,mouse_celltype, gene_labels ){

  human_celltype_stats = hu_markers |> filter(cell_type == human_celltype)
  mouse_celltype_stats = mou_markers |> filter(cell_type == mouse_celltype)

  gene_index = match(human_celltype_stats$gene, mouse_celltype_stats$gene)
  mouse_celltype_stats = mouse_celltype_stats[gene_index, ]

  #Highlight human specific genes
  # human_high = human_celltype_stats |> filter(log2(fold_change) >= 2) |> pull(gene)
  # mouse_low = mouse_celltype_stats |> filter(log2(fold_change) <= 0) |> pull(gene)
  # human_specific = intersect(human_high, mouse_low)

  #gene_labels_fc = c(gene_labels, human_specific)
  gene_labels_fc = gene_labels
  # human_high = human_celltype_stats |> filter(auroc <= .625 & auroc >= .375) |> pull(gene)
  # mouse_low = mouse_celltype_stats |> filter(auroc >= .75) |> pull(gene)
  # human_specific = intersect(human_high, mouse_low)

  #gene_labels_auroc = c(gene_labels, human_specific)

  gene_labels_auroc = gene_labels

  human_data = human_celltype_stats |> pull(stat_of_interest)
  mouse_data = mouse_celltype_stats |> pull(stat_of_interest)
  gene_label_vec_fc = human_celltype_stats$gene
  gene_label_vec_fc[!gene_label_vec_fc %in% gene_labels_fc] = NA
  color_vec_fc = rep('Not labeled', length = length(gene_label_vec_fc))
  color_vec_fc[!is.na(gene_label_vec_fc)] = 'Labeled'
  gene_label_vec_auroc = human_celltype_stats$gene
  gene_label_vec_auroc[!gene_label_vec_auroc %in% gene_labels_auroc] = NA
  color_vec_auroc = rep('Not labeled', length = length(gene_label_vec_auroc))
  color_vec_auroc[!is.na(gene_label_vec_auroc)] = 'Labeled'

  data_df = data.frame(human_celltype = human_celltype, mouse_celltype = mouse_celltype,
    human_stat = human_data, mouse_stat = mouse_data, 
    stat = stat_of_interest, gene_label_fc = gene_label_vec_fc, point_color_fc = color_vec_fc,
    gene_label_auroc = gene_label_vec_auroc, point_color_auroc = color_vec_auroc)

  if(stat_of_interest == 'auroc'){
    p = ggplot(data_df, aes(x = human_stat, y = mouse_stat, label = gene_label_auroc, 
      color = point_color_auroc, size = point_color_auroc, alpha = point_color_auroc)) + 
      # Mouse-specific box (x in [.375, .625], y >= .625)
      annotate("rect", xmin = 0.375, xmax = 0.625, ymin = 0.625, ymax = 1, 
               fill = "blue", alpha = 0.2) +
      # Human-specific box (y in [.375, .625], x >= .625)
      annotate("rect", xmin = 0.625, xmax = 1, ymin = 0.375, ymax = 0.625, 
               fill = "green", alpha = 0.2) +
      # Concordant box (x >= .625, y >= .625, excluding the other regions)
      annotate("rect", xmin = 0.625, xmax = 1, ymin = 0.625, ymax = 1, 
               fill = "purple", alpha = 0.2) +
      geom_point(alpha = .5, show.legend = FALSE) +
      geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') + 
      geom_label_repel(max.overlaps = Inf, min.segment.length = 0, size = 3) +
      scale_color_manual(values = c('Not labeled' = 'black', 'Labeled' = 'red')) +
      scale_size_manual(values = c('Not labeled' = 1, 'Labeled' = 2)) +
      scale_alpha_manual(values = c('Not labeled' = .25, 'Labeled' = 1) ) +
      coord_fixed() +
      xlab(sprintf('Human AUROC: %s', human_celltype)) + ylab(sprintf('Mouse AUROC: %s', mouse_celltype)) + ylim(0,1) + xlim(0,1) +
      theme_bw() + ggtitle(sprintf('%s', stat_of_interest))
    return(p)

  }else if (stat_of_interest == 'fold_change'){
    p = ggplot(data_df, aes(x = log2(human_stat), y = log2(mouse_stat), label = gene_label_fc,
    color = point_color_fc, size = point_color_fc, alpha = point_color_fc)) + 
      geom_point(alpha = .5, show.legened = FALSE) +
      geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') + 
      geom_label_repel(max.overlaps = Inf, min.segment.length = 0, size = 3) +
      scale_color_manual(values = c('Not labeled' = 'black', 'Labeled' = 'red')) +
      scale_size_manual(values = c('Not labeled' = 1, 'Labeled' = 2)) +
      scale_alpha_manual(values = c('Not labeled' = .25, 'Labeled' = 1) ) +
      xlab(sprintf('Human log2(FC): %s', human_celltype)) + ylab(sprintf('Mouse log2(FC): %s', mouse_celltype)) +
      theme_bw() + ggtitle(sprintf('%s', stat_of_interest))
    return(p)
  }
}


subP_auroc_p = plot_de_stats(human_markers, mouse_markers, 'auroc', 
human_celltype = 'MHb_A', mouse_celltype = 'MHb_A', gene_labels = c( 'SLC12A5', 'TAC1'))
subP_auroc_p
subP_fc_p = plot_de_stats(human_markers, mouse_markers, 'fold_change', 
human_celltype = 'MHb_A', mouse_celltype = 'MHb_A', gene_labels = c('SLC12A5', 'TAC1'))
subP_fc_p

chol_auroc_p = plot_de_stats(human_markers, mouse_markers, 'auroc', 
human_celltype = 'MHb_B', mouse_celltype = 'MHb_B', gene_labels = c('SLC12A5', 'SLC5A7'))
chol_auroc_p
chol_fc_p = plot_de_stats(human_markers, mouse_markers, 'fold_change', 
human_celltype = 'MHb_B', mouse_celltype = 'MHb_B', gene_labels = c('SLC12A5', 'SLC5A7'))
chol_fc_p


LHb_A_auroc_p = plot_de_stats(human_markers, mouse_markers, 'auroc', 
human_celltype = 'LHb_A', mouse_celltype = 'LHb_A', gene_labels = c('SLC12A5', 'HTR2C', 'CBLN2'))
LHb_A_auroc_p
LHb_A_fc_p = plot_de_stats(human_markers, mouse_markers, 'fold_change', 
human_celltype = 'LHb_A', mouse_celltype = 'LHb_A', gene_labels = c('SLC12A5', 'HTR2C', 'CBLN2'))
LHb_A_fc_p

LHb_B_auroc_p = plot_de_stats(human_markers, mouse_markers, 'auroc', 
human_celltype = 'LHb_B', mouse_celltype = 'LHb_B', gene_labels = c('SLC12A5', 'TENM1'))
LHb_B_auroc_p
LHb_B_fc_p = plot_de_stats(human_markers, mouse_markers, 'fold_change', 
human_celltype = 'LHb_B', mouse_celltype = 'LHb_B', gene_labels = c('SLC12A5', 'TENM1'))
LHb_B_fc_p

LHb_C_auroc_p = plot_de_stats(human_markers, mouse_markers, 'auroc', 
human_celltype = 'LHb_C', mouse_celltype = 'LHb_C', gene_labels = c('SLC12A5',  'GABRB1'))
LHb_C_auroc_p
LHb_C_fc_p = plot_de_stats(human_markers, mouse_markers, 'fold_change', 
human_celltype = 'LHb_C', mouse_celltype = 'LHb_C', gene_labels = c('SLC12A5','GABRB1'))
LHb_C_fc_p

table(human_markers$cell_type)



human_kcc2_stats = human_markers |>
  filter(gene == 'SLC12A5' & cell_type %in% c('MHb_A','MHb_B','LHb_A','LHb_B','LHb_C')) |> 
  select(cell_type, auroc) |> 
  mutate(species = 'Human', gene = 'SLC12A5')

mouse_kcc2_stats = mouse_markers |>
  filter(gene == 'SLC12A5' & cell_type %in% c('MHb_A','MHb_B','LHb_A','LHb_B','LHb_C')) |> 
  select(cell_type, auroc) |> 
  mutate(species = 'Mouse', gene = 'SLC12A5')

# Combine the data for grouped lollipop plot
plot_data <- bind_rows(human_kcc2_stats, mouse_kcc2_stats)
plot_data$cell_type = factor(plot_data$cell_type, levels = c('LHb_C','LHb_B','LHb_A','MHb_A','MHb_B'))
# Create grouped lollipop chart
lolly_plot = ggplot(plot_data, aes(x = cell_type, y = auroc, group = species)) +
  geom_line(aes(group = cell_type, color = cell_type), linewidth = 2) +
  scale_color_manual(values = my_colors_mid, name = "Cell Type") +
  ggnewscale::new_scale_color() +
  geom_point(size = 4, aes(color = species)) +
  scale_color_manual(values = c("Human" = "#e92f16", "Mouse" = "#22c3e7"), name = "Species") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  labs(
    x = "Cell Type",
    y = "AUROC",
    title = "KCC2 (SLC12A5) Marker DE by Cell Type"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
lolly_plot


#save plots

ggsave(plot = subP_auroc_p, filename = paste0(plot_path, '/subP_auroc_cross_species.pdf'), 
width = 6, height = 5, device = 'pdf')
ggsave(plot = chol_auroc_p, filename = paste0(plot_path, '/chol_auroc_cross_species.pdf'), 
width = 6, height = 5, device = 'pdf')
ggsave(plot = LHb_A_auroc_p, filename = paste0(plot_path, '/LHb_A_auroc_cross_species.pdf'), 
width = 6, height = 5, device = 'pdf')
ggsave(plot = LHb_B_auroc_p, filename = paste0(plot_path, '/LHb_B_auroc_cross_species.pdf'), 
width = 6, height = 5, device = 'pdf')
ggsave(plot = LHb_C_auroc_p, filename = paste0(plot_path, '/LHb_C_auroc_cross_species.pdf'), 
width = 6, height = 5, device = 'pdf')

ggsave(plot = lolly_plot, filename = paste0(plot_path, '/KCC2_lollipop_plot.pdf'), 
width = 6, height = 5, device = 'pdf')











