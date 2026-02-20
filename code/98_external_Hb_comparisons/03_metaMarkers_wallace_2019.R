#Script for getting MetaMarkers from the metaClusters across individual mouse donors from the Wallace 2019 dataset
#Will likely be updating metaCluster annotations, just setting up the code to check out the markers for those clusters


library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
library(dplyr)
library(ggplot2)
library(here)

here::here()

#Path to the Wallace 2019 data
wallace_path = here('processed-data', '98_external_Hb_comparisons', 'Wallace_etal_2019_habenula_scseq')
list.files(wallace_path)
#Path to save any generated data
new_data_path = here('processed-data', '98_external_Hb_comparisons', '03_metaMarkers_wallace_2019')
#Path to already generated data
prev_data_path = here('processed-data', '98_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')
#Path to plot directory
plot_path = here('plots', '98_external_Hb_comparisons', '03_metaMarkers_wallace_2019')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)

#Load up the full SCE object, contains the metacluster annotations
all_donor_sce = readRDS(paste0(prev_data_path, '/all_donor_sce_with_denovo_clusters.rds'))

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
export_markers(markers_160822, paste0(new_data_path, '/markers_160822_meta_clusters_markers.csv'))
export_markers(markers_161103, paste0(new_data_path, '/markers_161103_meta_clusters_markers.csv'))
export_markers(markers_161102, paste0(new_data_path, '/markers_161102_meta_clusters_markers.csv'))
export_markers(markers_161105, paste0(new_data_path, '/markers_161105_meta_clusters_markers.csv'))


#Load up markers and get the metaMarkers 
wallace_mouse_hab_markers = list(
    donor_160822 = read_markers(paste0(new_data_path, '/markers_160822_meta_clusters_markers.csv.gz')),
    donor_161102 = read_markers(paste0(new_data_path, '/markers_161102_meta_clusters_markers.csv.gz')),
    donor_161103 = read_markers(paste0(new_data_path, '/markers_161103_meta_clusters_markers.csv.gz')),
    donor_161105 = read_markers(paste0(new_data_path, '/markers_161105_meta_clusters_markers.csv.gz'))
    
)


wallace_mouse_hab_metaM = make_meta_markers(wallace_mouse_hab_markers, detailed_stats = TRUE)

#Save the metamarkers
export_meta_markers(wallace_mouse_hab_metaM, 
  paste0(new_data_path, '/mouse_hab_meta_markers.csv'), 
  names(wallace_mouse_hab_metaM))

wallace_mouse_hab_metaM = read_meta_markers(paste0(new_data_path, '/mouse_hab_meta_markers.csv.gz'))

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
indv_mouse_seurats = readRDS(file = paste0(prev_data_path, '/individual_donor_seurat_objects_list.rds'))

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


#Markers used in the original Wallace 2019 paper in Figure 1
custom_markers = c('Tac2', 'Slc17a7', 'Slc17a6', 'Snap25', 'Gap43', 'Slc6a11', 'Cldn5', 'Abcc9', 'Pdgfrb', 'Cx3cr1', 'Mrc1', 'Col3a1', 'Gpr17', 'Mog', 'Olig1', 'Pdgfra')


p_bubble = get_bubble_plot(hab_160822_seurat, custom_markers, 'Mouse: 160822')
p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(hab_161102_seurat, custom_markers, 'Mouse: 161102')
p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(hab_161103_seurat, custom_markers, 'Mouse: 161103')
p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(hab_161105_seurat, custom_markers, 'Mouse: 161105')
p_bubble[[1]]
p_bubble[[2]]

#Human marker panel
#We use Tac3 in humans, but not present in the mouse annotations
custom_markers = c('Tac2','Gpr151','Pou4f1','Mbp')


p_bubble = get_bubble_plot(hab_160822_seurat, custom_markers, 'Mouse: 160822')
p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(hab_161102_seurat, custom_markers, 'Mouse: 161102')
p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(hab_161103_seurat, custom_markers, 'Mouse: 161103')
p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(hab_161105_seurat, custom_markers, 'Mouse: 161105')
p_bubble[[1]]
p_bubble[[2]]


#Neuronal subtype markers
#Chat, Slc18a3, Slc5a7 are all cholinergic markers
#Tac1 is a marker for substance P neurons
#Slc17a7 and a6 are Vglut1 and Vglut2 
custom_markers = c('Chat', 'Slc18a3', 'Slc5a7','Tac1', 'Slc17a7', 'Slc17a6')


p_bubble = get_bubble_plot(hab_160822_seurat, custom_markers, 'Mouse: 160822')
p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(hab_161102_seurat, custom_markers, 'Mouse: 161102')
p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(hab_161103_seurat, custom_markers, 'Mouse: 161103')
p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(hab_161105_seurat, custom_markers, 'Mouse: 161105')
p_bubble[[1]]
p_bubble[[2]]


#And repeat the above but looking at the author habenula subclusters
p_bubble = get_bubble_plot(hab_160822_seurat, custom_markers, 'Mouse: 160822', group_col = 'author_subHab_celltype')
p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(hab_161102_seurat, custom_markers, 'Mouse: 161102', group_col = 'author_subHab_celltype')
p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(hab_161103_seurat, custom_markers, 'Mouse: 161103', group_col = 'author_subHab_celltype')
p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(hab_161105_seurat, custom_markers, 'Mouse: 161105', group_col = 'author_subHab_celltype')
p_bubble[[1]]
p_bubble[[2]]



#I think we can confidently label the non-neuronal metaclusters, the cholinergic - supstance P medial Hab clusters, and at least the Lateral hab clusters
#Save bubble plots using the full data across all donors
#And save the metadata with the annotated MetaClusters, use this moving forward

#Switch the full dataset to Seurat for the bubble plots
all_donor_seurat = as.Seurat(all_donor_sce, counts = "counts", data = "cpm")
all_donor_seurat

#Curious to see the meta-cluster annotations in the full dataset
#repeat the standard dim reduction from seurat

#HVGs
all_donor_seurat <- FindVariableFeatures(all_donor_seurat, selection.method = "vst", nfeatures = 2000)

#Scale data
all.genes <- rownames(all_donor_seurat)
all_donor_seurat <- ScaleData(all_donor_seurat, features = all.genes)

#PCA
all_donor_seurat  <- RunPCA(all_donor_seurat , features = VariableFeatures(object = all_donor_seurat ))
DimPlot(all_donor_seurat, reduction = "pca") + NoLegend()

#UMAP
all_donor_seurat  <- RunUMAP(all_donor_seurat , dims = 1:20)
p1 = DimPlot(all_donor_seurat , reduction = "umap", group.by = 'meta_cluster', label = TRUE) + 
  ggtitle('Wallace 2019: MetaCluster annotations')
p1
ggsave(p1, filename = 'wallace_mouse_hab_meta_cluster_umap.pdf', path = plot_path,
device = 'pdf', width = 8, height = 7)

#While we're here, save donor umap, author annotation umap
p2 = DimPlot(all_donor_seurat , reduction = "umap", group.by = 'putative_donor') + 
  ggtitle('Wallace 2019: mouse sample annotations')
p2
ggsave(p2, filename = 'wallace_mouse_hab_donor_umap.pdf', path = plot_path,
device = 'pdf', width = 8, height = 7)


p3 = DimPlot(all_donor_seurat , reduction = "umap", group.by = 'author_celltype', label = TRUE) + 
  ggtitle('Wallace 2019: author annotations')
p3
ggsave(p3, filename = 'wallace_mouse_hab_author_annot_umap.pdf', path = plot_path,
device = 'pdf', width = 8, height = 7)

p4 = DimPlot(all_donor_seurat , reduction = "umap", group.by = 'author_subHab_celltype', label = TRUE) + 
  ggtitle('Wallace 2019: author annotations')
p4
ggsave(p4, filename = 'wallace_mouse_hab_author_subHab_annot_umap.pdf', path = plot_path,
device = 'pdf', width = 8, height = 7)


#And now the bubble plot with the Wallace marker panel
custom_markers = c('Tac2', 'Slc17a7', 'Slc17a6', 'Snap25', 'Gap43', 'Slc6a11', 'Cldn5', 'Abcc9', 'Pdgfrb', 'Cx3cr1', 'Mrc1', 'Col3a1', 'Gpr17', 'Mog', 'Olig1', 'Pdgfra')

p_bubble_paper_markers = get_bubble_plot(all_donor_seurat , custom_markers, 'Wallace 2019: all mouse samples')
p_bubble_paper_markers
ggsave(p_bubble_paper_markers[[2]], filename = 'wallace_mouse_hab_paper_markers_meta_cluster_bubble.pdf', path = plot_path,
device = 'pdf', width = 10, height = 8)


#And the bubble plot with the human marker panel
custom_markers = c('Tac2','Gpr151','Pou4f1','Mbp')

p_bubble_human_markers = get_bubble_plot(all_donor_seurat , custom_markers, 'Wallace 2019: all mouse samples')
p_bubble_human_markers
ggsave(p_bubble_human_markers[[2]], filename = 'wallace_mouse_human_hab_markers_meta_cluster_bubble.pdf', path = plot_path,
device = 'pdf', width = 10, height = 8)


#And the cholinergic and substance P markers
custom_markers = c('Chat', 'Slc18a3', 'Slc5a7','Tac1', 'Slc17a7', 'Slc17a6')

p_bubble_excit_markers = get_bubble_plot(all_donor_seurat , custom_markers, 'Wallace 2019: all mouse samples')
p_bubble_excit_markers
ggsave(p_bubble_excit_markers[[2]], filename = 'wallace_mouse_excite_subtype_markers_meta_cluster_bubble.pdf', path = plot_path,
device = 'pdf', width = 10, height = 8)


#metacluster annotations from all the above
meta_annot_vec = c('MHb_subP' = 'meta_cluster1',
                    'Polydendrocytes' = 'meta_cluster2',
                    'Endothelial' = 'meta_cluster3',
                    'LHb_1' = 'meta_cluster4',
                    'Macrophages' = 'meta_cluster5',
                    'Astrocytes' = 'meta_cluster6',
                    'Fibroblasts' = 'meta_cluster7',
                    'MHb_cholinergic' = 'meta_cluster8',
                    'Oligodendrocytes' = 'meta_cluster9',
                    'MHb_subP_cholinergic' = 'meta_cluster10',
                    'Differentiating Oligodendrocytes' = 'meta_cluster11',
                    'Pericytes' = 'meta_cluster12',
                    'Microglia' = 'meta_cluster13',
                    'LHb_2' = 'meta_cluster14',
                    'outliers' = 'outliers'  
)

meta_annot_vec  = setNames(names(meta_annot_vec), meta_annot_vec)


all_donor_seurat$meta_clust_celltype_annot = unname(meta_annot_vec[all_donor_seurat$meta_cluster])
table(all_donor_seurat$meta_clust_celltype_annot, all_donor_seurat$meta_cluster)

#Save the metadata as a data.frame to add to the seurat data object later
full_seurat_metadata = all_donor_seurat@meta.data
saveRDS(full_seurat_metadata, paste0(new_data_path, '/wallace_mouse_metaclust_celltype_annot_metadata.rds'))



#and a final umap with the annotated metaclusters
p5 = DimPlot(all_donor_seurat , reduction = "umap", group.by = 'meta_clust_celltype_annot', label = TRUE) + 
  ggtitle('Wallace 2019: metacluster celltype annotations')
p5
ggsave(p5, filename = 'wallace_mouse_hab_metaCluster_celltype_annot_umap.pdf', path = plot_path,
device = 'pdf', width = 8, height = 7)








