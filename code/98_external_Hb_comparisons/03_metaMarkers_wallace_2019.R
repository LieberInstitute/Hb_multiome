#Script for getting MetaMarkers from the metaClusters across individual mouse donors from the Wallace 2019 dataset
#Will likely be updating metaCluster annotations, just setting up the code to check out the markers for those clusters


library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
library(dplyr)
library(ggplot2)


#Load up the full SCE object, contains the metacluster annotations
all_donor_sce = readRDS('processed-data/98_external_Hb_comparisons/02_qc_and_clust_wallace_2019/all_donor_sce_with_denovo_clusters.rds')


#Generate MetaMarkers from the metaclusters
#Split into individual donor SCE objects and compute intial markers for the metaclusters
table(all_donor_sce$putative_donor)

sce_160822 <- all_donor_sce[, all_donor_sce$putative_donor == "160822"]
sce_161102 <- all_donor_sce[, all_donor_sce$putative_donor == "161102"]
sce_161103 <- all_donor_sce[, all_donor_sce$putative_donor == "161103"]
sce_161105 <- all_donor_sce[, all_donor_sce$putative_donor == "161105"]

markers_160822 = compute_markers(assay(sce_160822, "cpm"), sce_160822$meta_cluster)
markers_161102 = compute_markers(assay(sce_161102, "cpm"), sce_161102$meta_cluster)
markers_161103 = compute_markers(assay(sce_161103, "cpm"), sce_161103$meta_cluster)
markers_161105 = compute_markers(assay(sce_161105, "cpm"), sce_161105$meta_cluster)

head(markers_160822)

#Save markers
export_markers(markers_160822, 'processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_160822_meta_clusters_markers.csv')
export_markers(markers_161103, 'processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_161103_meta_clusters_markers.csv')
export_markers(markers_161102, 'processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_161102_meta_clusters_markers.csv')
export_markers(markers_161105, 'processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_161105_meta_clusters_markers.csv')


#Load up markers and get the metaMarkers 
wallace_mouse_hab_markers = list(
    donor_160822 = read_markers("processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_160822_meta_clusters_markers.csv.gz"),
    donor_161102 = read_markers("processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_161102_meta_clusters_markers.csv.gz"),
    donor_161103 = read_markers("processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_161103_meta_clusters_markers.csv.gz"),
    donor_161105 = read_markers("processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_161105_meta_clusters_markers.csv.gz")
    
)


wallace_mouse_hab_metaM = make_meta_markers(wallace_mouse_hab_markers, detailed_stats = TRUE)

#Save the metamarkers
export_meta_markers(wallace_mouse_hab_metaM, 
  "processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/mouse_hab_meta_markers.csv", 
  names(wallace_mouse_hab_metaM))

wallace_mouse_hab_metaM = read_meta_markers("processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/mouse_hab_meta_markers.csv.gz")

wallace_mouse_hab_metaM %>% group_by(cell_type) %>% slice_min(rank, n = 20) %>% View()

wallace_mouse_hab_metaM %>% filter(gene %in% c('Chat','Slc18a3','Slc5a7', 'Tac1', 'Vglut1', 'Vglut2')) %>% View()





#Check out some of the top markers used in the Wallace 2019 paper
#Tac2, Slc17a7 are MHb markers
#Slc17a6, Snap25 are MHb and LHb markers
#Gap43 is an LHb marker
wallace_mouse_hab_metaM %>% filter(gene %in% c('Tac2', 'Slc17a7', 'Slc17a6', 'Snap25', 'Gap43')) %>%
  group_by(cell_type) %>% arrange(rank, .by_group = T) %>% View()



#Astrocyte, endothelial, pericyte, pericyte
wallace_mouse_hab_metaM %>% filter(gene %in% c('Slc6a11', 'Cldn5', 'Abcc9', 'Pdgfrb')) %>%
  group_by(cell_type) %>% arrange(rank, .by_group = T) %>% View()



#Microglia, macrophages, fibroblasts
wallace_mouse_hab_metaM %>% filter(gene %in% c('Cx3cr1', 'Mrc1', 'Col3a1')) %>%
  group_by(cell_type) %>% arrange(rank, .by_group = T) %>% View()


#Mog is an oligo and Diff. Oligo marker 
#GPR17 is a Diff. Oligo and Polydendrocyte marker
#Olig1 is a marker for all three
#Pdgfra is a marker for polydendrocytes
wallace_mouse_hab_metaM %>% filter(gene %in% c('Gpr17', 'Mog', 'Olig1', 'Pdgfra')) %>%
  group_by(cell_type) %>% arrange(rank, .by_group = T) %>% View()




#Load up the individual seurat objects, add the metaCluster annotations, and check out what it looks like in the UMAPs
indv_mouse_seurats = readRDS(file = 'processed-data/98_external_Hb_comparisons/02_qc_and_clust_wallace_2019/individual_donor_seurat_objects_list.rds')

hab_160822_seurat = indv_mouse_seurats$hab_160822
hab_161102_seurat = indv_mouse_seurats$hab_161102
hab_161103_seurat = indv_mouse_seurats$hab_161103
hab_161105_seurat = indv_mouse_seurats$hab_161105


full_cell_names_160822 = paste('hab_160822', hab_160822_seurat$seurat_clusters, sep = '|')
hab_160822_seurat$meta_cluster = unname(all_donor_sce$meta_cluster[full_cell_names_160822])

full_cell_names_161102 = paste('hab_161102', hab_161102_seurat$seurat_clusters, sep = '|')
hab_161102_seurat$meta_cluster = unname(all_donor_sce$meta_cluster[full_cell_names_161102])

full_cell_names_161103 = paste('hab_161103', hab_161103_seurat$seurat_clusters, sep = '|')
hab_161103_seurat$meta_cluster = unname(all_donor_sce$meta_cluster[full_cell_names_161103])

full_cell_names_161105 = paste('hab_161105', hab_161105_seurat$seurat_clusters, sep = '|')
hab_161105_seurat$meta_cluster = unname(all_donor_sce$meta_cluster[full_cell_names_161105])



DimPlot(hab_160822_seurat, group.by = 'meta_cluster', label = TRUE) + ggtitle('Mouse: 160822')
DimPlot(hab_161102_seurat, group.by = 'meta_cluster', label = TRUE) + ggtitle('Mouse: 161102')
DimPlot(hab_161103_seurat, group.by = 'meta_cluster', label = TRUE) + ggtitle('Mouse: 161103')
DimPlot(hab_161105_seurat, group.by = 'meta_cluster', label = TRUE) + ggtitle('Mouse: 161105')



#Check out some bubble plots of the top metaMarkers across the metaclusters
top_cluster_markers = wallace_mouse_hab_metaM %>% 
  filter(rank <= 5) %>% pull(gene)
top_cluster_markers = unique(top_cluster_markers)

DotPlot(hab_160822_seurat, 
        features = c(top_cluster_markers),
        group.by = 'meta_cluster')

DotPlot(hab_161105_seurat, 
        features = c(top_cluster_markers),
        group.by = 'meta_cluster')



#Custom bubble plot function

get_bubble_plot = function(seurat_object, top_markers, sample_name, group_col = 'meta_cluster'){
  # Extract expression data and metadata
  expr_data <- FetchData(seurat_object, vars = top_markers, slot = "data")
  metadata <- seurat_object@meta.data

  # Combine into a data frame
  plot_data <- cbind(expr_data, meta_cluster = metadata[[group_col]]) %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = -meta_cluster, names_to = "gene", values_to = "expression")

  # Calculate mean expression and percent expressing per cluster
  summary_data <- plot_data %>%
    group_by(gene, meta_cluster) %>%
    summarise(
      mean_expression = mean(expression),
      pct_expressing = sum(expression > 0) / n() * 100,
      .groups = "drop"
    )

  # Set factor levels to control axis order
  summary_data$gene <- factor(summary_data$gene, levels = top_markers)
  summary_data$meta_cluster <- factor(summary_data$meta_cluster, 
                                      levels = sort(unique(summary_data$meta_cluster)))

  # Create bubble plot
  p1 = ggplot(summary_data, aes(x = gene, y = meta_cluster, size = mean_expression, color = pct_expressing)) +
    geom_point() +
    scale_color_gradient(low = "lightgrey", high = "red") +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() + ggtitle(sample_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = "Meta Cluster", size = "Mean Expression", color = "% Expressing")
  return(p1)
}


#Markers used in the original Wallace 2019 paper in Figure 1
custom_markers = c('Tac2', 'Slc17a7', 'Slc17a6', 'Snap25', 'Gap43', 'Slc6a11', 'Cldn5', 'Abcc9', 'Pdgfrb', 'Cx3cr1', 'Mrc1', 'Col3a1', 'Gpr17', 'Mog', 'Olig1', 'Pdgfra')


p_bubble = get_bubble_plot(hab_160822_seurat, custom_markers, 'Mouse: 160822')
p_bubble

p_bubble = get_bubble_plot(hab_161102_seurat, custom_markers, 'Mouse: 161102')
p_bubble

p_bubble = get_bubble_plot(hab_161103_seurat, custom_markers, 'Mouse: 161103')
p_bubble

p_bubble = get_bubble_plot(hab_161105_seurat, custom_markers, 'Mouse: 161105')
p_bubble


#Human marker panel
#We use Tac3 in humans, but not present in the mouse annotations
custom_markers = c('Tac2','Gpr151','Pou4f1','Mbp')


p_bubble = get_bubble_plot(hab_160822_seurat, custom_markers, 'Mouse: 160822')
p_bubble

p_bubble = get_bubble_plot(hab_161102_seurat, custom_markers, 'Mouse: 161102')
p_bubble

p_bubble = get_bubble_plot(hab_161103_seurat, custom_markers, 'Mouse: 161103')
p_bubble

p_bubble = get_bubble_plot(hab_161105_seurat, custom_markers, 'Mouse: 161105')
p_bubble


#Neuronal subtype markers
#Chat, Slc18a3, Slc5a7 are all cholinergic markers
#Tac1 is a marker for substance P neurons
#Slc17a7 and a6 are Vglut1 and Vglut2 
custom_markers = c('Chat', 'Slc18a3', 'Slc5a7','Tac1', 'Slc17a7', 'Slc17a6')


p_bubble = get_bubble_plot(hab_160822_seurat, custom_markers, 'Mouse: 160822')
p_bubble

p_bubble = get_bubble_plot(hab_161102_seurat, custom_markers, 'Mouse: 161102')
p_bubble

p_bubble = get_bubble_plot(hab_161103_seurat, custom_markers, 'Mouse: 161103')
p_bubble

p_bubble = get_bubble_plot(hab_161105_seurat, custom_markers, 'Mouse: 161105')
p_bubble


#And repeat the above but looking at the author habenula subclusters
p_bubble = get_bubble_plot(hab_160822_seurat, custom_markers, 'Mouse: 160822', group_col = 'author_subHab_celltype')
p_bubble

p_bubble = get_bubble_plot(hab_161102_seurat, custom_markers, 'Mouse: 161102', group_col = 'author_subHab_celltype')
p_bubble

p_bubble = get_bubble_plot(hab_161103_seurat, custom_markers, 'Mouse: 161103', group_col = 'author_subHab_celltype')
p_bubble

p_bubble = get_bubble_plot(hab_161105_seurat, custom_markers, 'Mouse: 161105', group_col = 'author_subHab_celltype')
p_bubble






