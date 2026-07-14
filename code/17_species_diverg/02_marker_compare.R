#Comparing marker stats between human and mouse for the consensus cell-type annotations

library(ggplot2)
library(ggrepel)
library(dplyr)
library(here)
library(MetaMarkers)

here::here()


#Path to save any generated data
new_data_path = here('processed-data', '05_03_annotation_adjustments', '15_marker_GO_enrich')
#Path to plot directory
plot_path = here('plots','05_03_annotation_adjustments', '15_marker_GO_enrich')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


human_marker_path = here('processed-data', '05_03_annotation_adjustments', '08_metamarkers')
mouse_marker_path = here('processed-data','17_species_diverg','01_hashikawa_markers')


#human_marker_stats = readRDS(file = paste0(human_marker_path, '/marker_stats_combo.rds'))
#mouse_marker_stats = readRDS(file = paste0(mouse_marker_path, '/marker_stats_combo.rds'))

human_markers = read_meta_markers(paste0(human_marker_path, '/multiome_refined_mid_meta_markers.csv.gz'))
mouse_markers = read_meta_markers(paste0(mouse_marker_path, '/hashikawa_meta_markers.csv.gz'))

#Filter for shared genes
shared_genes = intersect(unique(human_markers$gene), unique(mouse_markers$gene))
length(shared_genes)

human_markers = human_markers[human_markers$gene %in% shared_genes, ]
mouse_markers = mouse_markers[mouse_markers$gene %in% shared_genes, ]


human_markers
mouse_markers


plot_de_stats = function(hu_markers, mou_markers, stat_of_interest, celltype, gene_labels ){

  human_celltype_stats = hu_markers |> filter(cell_type == celltype)
  mouse_celltype_stats = mou_markers |> filter(cell_type == celltype)

  gene_index = match(human_celltype_stats$gene, mouse_celltype_stats$gene)
  mouse_celltype_stats = mouse_celltype_stats[gene_index, ]

  #Highlight human specific genes
  human_high = human_celltype_stats |> filter(log2(fold_change) >= 2) |> pull(gene)
  mouse_low = mouse_celltype_stats |> filter(log2(fold_change) <= 0) |> pull(gene)
  human_specific = intersect(human_high, mouse_low)

  gene_labels_fc = c(gene_labels, human_specific)
  
  human_high = human_celltype_stats |> filter(auroc >= .75) |> pull(gene)
  mouse_low = mouse_celltype_stats |> filter(auroc < .75) |> pull(gene)
  human_specific = intersect(human_high, mouse_low)

  gene_labels_auroc = c(gene_labels, human_specific)


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

  data_df = data.frame(celltype = celltype, human_stat = human_data, mouse_stat = mouse_data, 
    stat = stat_of_interest, gene_label_fc = gene_label_vec_fc, point_color_fc = color_vec_fc,
    gene_label_auroc = gene_label_vec_auroc, point_color_auroc = color_vec_auroc)

  if(stat_of_interest == 'auroc'){
    p = ggplot(data_df, aes(x = human_stat, y = mouse_stat, label = gene_label_auroc, 
      color = point_color_auroc, size = point_color_auroc, alpha = point_color_auroc)) + 
      geom_point(alpha = .5, show.legend = FALSE) +
      geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') + 
      geom_label_repel(max.overlaps = Inf, min.segment.length = 0) +
      scale_color_manual(values = c('Not labeled' = 'black', 'Labeled' = 'red')) +
      scale_size_manual(values = c('Not labeled' = 1, 'Labeled' = 2)) +
      scale_alpha_manual(values = c('Not labeled' = .25, 'Labeled' = 1) ) +
      coord_fixed() +
      xlab('Human AUROC') + ylab('Mouse AUROC') + ylim(0,1) + xlim(0,1) +
      theme_bw() + ggtitle(sprintf('%s: %s',celltype, stat_of_interest))
    return(p)

  }else if (stat_of_interest == 'fold_change'){
    p = ggplot(data_df, aes(x = log2(human_stat), y = log2(mouse_stat), label = gene_label_fc,
    color = point_color_fc, size = point_color_fc, alpha = point_color_fc)) + 
      geom_point(alpha = .5, show.legened = FALSE) +
      geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') + 
      geom_label_repel(max.overlaps = Inf, min.segment.length = 0, size = 5) +
      scale_color_manual(values = c('Not labeled' = 'black', 'Labeled' = 'red')) +
      scale_size_manual(values = c('Not labeled' = 1, 'Labeled' = 2)) +
      scale_alpha_manual(values = c('Not labeled' = .25, 'Labeled' = 1) ) +
      xlab('Human log2(FC)') + ylab('Mouse log2(FC)') +
      theme_bw() + ggtitle(sprintf('%s: %s',celltype, stat_of_interest))
    return(p)
  }
}

astro_auroc_p = plot_de_stats(human_markers, mouse_markers, 'auroc', 'Astrocyte', gene_labels = c('AQP4'))
astro_auroc_p
astro_fc_p = plot_de_stats(human_markers, mouse_markers, 'fold_change', 'Astrocyte', gene_labels = c('AQP4'))
astro_fc_p


subP_auroc_p = plot_de_stats(human_markers, mouse_markers, 'auroc', 'MHb.1', gene_labels = c('OPRM1'))
subP_auroc_p
subP_fc_p = plot_de_stats(human_markers, mouse_markers, 'fold_change', 'MHb.1', gene_labels = 'OPRM1')
subP_fc_p

chol_auroc_p = plot_de_stats(human_markers, mouse_markers, 'auroc', 'MHb.2', gene_labels = c('CHRNA3','CHRNB4','CHRNA5'))
chol_auroc_p
chol_fc_p = plot_de_stats(human_markers, mouse_markers, 'fold_change', 'MHb.2', gene_labels = c('CHRNA3','CHRNB4','CHRNA5'))
chol_fc_p


LHb_A_auroc_p = plot_de_stats(human_markers, mouse_markers, 'auroc', 'LHb.2.7', gene_labels = c('OPRM1'))
LHb_A_auroc_p
LHb_A_fc_p = plot_de_stats(human_markers, mouse_markers, 'fold_change', 'LHb.2.7', gene_labels = c('OPRM1'))
LHb_A_fc_p

LHb_B_auroc_p = plot_de_stats(human_markers, mouse_markers, 'auroc', 'LHb.1.3.4', gene_labels = c('OPRM1'))
LHb_B_auroc_p
LHb_B_fc_p = plot_de_stats(human_markers, mouse_markers, 'fold_change', 'LHb.1.3.4', gene_labels = c('OPRM1'))
LHb_B_fc_p

LHb_C_auroc_p = plot_de_stats(human_markers, mouse_markers, 'auroc', 'LHb.4', gene_labels = c('OPRM1'))
LHb_C_auroc_p
LHb_C_fc_p = plot_de_stats(human_markers, mouse_markers, 'fold_change', 'LHb.4', gene_labels = c('OPRM1'))
LHb_C_fc_p

table(human_markers$cell_type)




















