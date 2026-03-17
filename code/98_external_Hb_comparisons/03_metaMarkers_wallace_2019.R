#Script for getting MetaMarkers from the metaClusters across individual mouse donors from the Wallace 2019 dataset
#Will likely be updating metaCluster annotations, just setting up the code to check out the markers for those clusters


library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
library(MetaNeighbor)
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


#Source the bubble plot functions
source(here('code','98_external_Hb_comparisons', 'bubble_plot_functions.R'))


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


#GLS and GLUD1/2 expression
custom_markers = c('Gls','Gls2', 'Glud1', 'Slc17a7', 'Slc17a6', 'Gad1','Gad2','Slc32a1', 'Pvalb')

p_bubble_excit_markers = get_bubble_plot(all_donor_seurat , custom_markers, 'Wallace 2019: all mouse samples', 
group_col = 'meta_clust_celltype_annot')
p_bubble_excit_markers


##########################
#Look for potential lateral habenula GABAergic cells
#########################


custom_markers = c('Gad1', 'Gad2', 'Slc32a1','Pvalb', 'Slc17a6', 'Slc17a7')

p_bubble_gaba_glut_markers = get_bubble_plot(all_donor_seurat , custom_markers, 'Wallace 2019: all mouse samples', 
group_col = 'meta_clust_celltype_annot')
p_bubble_gaba_glut_markers


#The LHb1 has the highest Gad2 expression and VGAT expression, but both are low in general and small percentages. This matches the trend in humans
mouse_gad1_p = FeaturePlot(all_donor_seurat, features = "Gad1", reduction = "umap", pt.size = 1, slot = 'scale.data') +
  scale_color_gradient(low = "white", high = "red", name = 'z-score')

mouse_gad2_p = FeaturePlot(all_donor_seurat, features = "Gad2", reduction = "umap", pt.size = 1, slot = 'scale.data') +
  scale_color_gradient(low = "white", high = "red", name = 'z-score')

mouse_pvalb_p = FeaturePlot(all_donor_seurat, features = "Pvalb", reduction = "umap", pt.size = 1, slot = 'scale.data') +
  scale_color_gradient(low = "white", high = "red", name = 'z-score')

mouse_vgat_p = FeaturePlot(all_donor_seurat, features = "Slc32a1", reduction = "umap", pt.size = 1, slot = 'scale.data') +
  scale_color_gradient(low = "white", high = "red", name = 'z-score')


#Look at the co-expression of specific genes
# Get expression data
umap_data <- as.data.frame(Embeddings(all_donor_seurat, reduction = "umap"))
umap_data$GAD2 <- FetchData(all_donor_seurat, vars = "Gad2", slot = "data")[, 1]
umap_data$PVALB <- FetchData(all_donor_seurat, vars = "Pvalb", slot = "data")[, 1]
umap_data$SLC17A6 <- FetchData(all_donor_seurat, vars = "Slc17a6", slot = "data")[, 1]

# Create a coexpression category
umap_data$coexpression <- ifelse(umap_data$GAD2 > 0 & umap_data$SLC17A6 > 0, "Both",
                                  ifelse(umap_data$GAD2 > 0, "Gad2 only",
                                         ifelse(umap_data$SLC17A6 > 0, "Slc17a6 only", "Neither")))

gad2_vglut2_p = ggplot(umap_data, aes(x = umap_1, y = umap_2, color = coexpression)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Both" = "purple", "Gad2 only" = "red", "Slc17a6 only" = "blue", "Neither" = "lightgrey")) +
  theme_bw() +
  labs(title = "Gad2 and Slc17a6 Co-expression: Wallace mouse")


umap_data$coexpression <- ifelse(umap_data$GAD2 > 0 & umap_data$PVALB > 0, "Both",
                                  ifelse(umap_data$GAD2 > 0, "Gad2 only",
                                         ifelse(umap_data$PVALB > 0, "Pvalb only", "Neither")))

gad2_pvalb_p = ggplot(umap_data, aes(x = umap_1, y = umap_2, color = coexpression)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Both" = "purple", "Gad2 only" = "red", "Pvalb only" = "blue", "Neither" = "lightgrey")) +
  theme_bw() +
  labs(title = "Gad2 and Pvalb Co-expression: Wallace mouse")


p_bubble_gaba_glut_markers[[1]]
p_bubble_gaba_glut_markers[[2]]

mouse_gad1_p 
mouse_gad2_p
mouse_pvalb_p
mouse_vgat_p

gad2_vglut2_p
gad2_pvalb_p 

ggsave(p_bubble_gaba_glut_markers[[1]], filename = 'wallace_mouse_gaba_glut_meta_annots_bubble.pdf', path = plot_path,
device = 'pdf', width = 10, height = 8)

ggsave(p_bubble_gaba_glut_markers[[2]], filename = 'wallace_mouse_gaba_glut_zscore_meta_annots_bubble.pdf', path = plot_path,
device = 'pdf', width = 10, height = 8)

ggsave(mouse_gad1_p , filename = 'wallace_mouse_gad1_exp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(mouse_gad2_p , filename = 'wallace_mouse_gad2_exp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(mouse_pvalb_p , filename = 'wallace_mouse_pvalb_exp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(mouse_vgat_p , filename = 'wallace_mouse_vgat_exp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(gad2_vglut2_p , filename = 'wallace_mouse_gad2_vglut2_coexp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(gad2_pvalb_p , filename = 'wallace_mouse_gad2_pvalb_coexp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

#Coexpression plots of the GABA-Glut genes

all_donor_sce$meta_clust_celltype_annot = unname(meta_annot_vec[all_donor_sce$meta_cluster])
logcounts(all_donor_sce) = log1p(assay(all_donor_sce, "cpm"))
#Plotting co-expression of excitatory and inhibitory markers

pairwise_coexpression <- function(mat, genes, cluster_name) {
  # mat: genes x cells matrix for one cluster
  
  detected <- mat[genes, , drop = FALSE] > 0
  
  res <- expand.grid(gene1 = genes, gene2 = genes, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(percent = mean(detected[gene1, ] & detected[gene2, ]) * 100) %>%
    ungroup() %>%
    mutate(cluster = cluster_name)
  
  res
}


#Full dataset. Of note, this does not use the corrected counts, the corrected counts were not saved with this version of the data
genes <- c('Slc32a1',"Gad1","Gad2","Slc17a6", "Slc17a7")
expr_mat <- assay(all_donor_sce, "logcounts")

cluster_to_annotate = "meta_clust_celltype_annot"
clusters <- unique(colData(all_donor_sce)[[cluster_to_annotate]])

coexp_df <- lapply(clusters, function(cl) {
  cells <- colData(all_donor_sce)[[cluster_to_annotate]] == cl
  mat_sub <- expr_mat[, cells, drop = FALSE]
  pairwise_coexpression(mat_sub, genes, cluster_name = cl)
}) %>%
  bind_rows()

coexp_df$gene1 <- factor(coexp_df$gene1, levels = genes)
coexp_df$gene2 <- factor(coexp_df$gene2, levels = rev(genes))


#Edited the original code to have grey be between 0-5%, previously the very low percentages were difficult to see.
p <- ggplot(coexp_df, aes(x = gene1, y = gene2, fill = percent)) +
  geom_tile(color = "grey70", linewidth = 0.3) +
  facet_wrap(~ cluster, nrow = 5) +
  scale_fill_gradientn(
    colours = c("grey95","#f1e2c6", "#f1e2c6", "#f0c94a", "#df8b27", "#d92523", "#8b0d19"),
    values = c(0, 0.05, 0.20, 0.40, 0.60, 0.8, 1),
    limits = c(0, 100),
    breaks = c( 5, 20, 40, 60, 80, 100),
    name = "Percent of cells expressing\ntwo genes"
  ) +
  coord_equal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  xlab(NULL) +
  ylab(NULL) + ggtitle("Co-expression of GABA and Glut markers in Wallace mouse dataset")

print(p)

ggsave(path = plot_path, filename = 'GABA_Glut_coexpression_wallace_mouse.pdf', plot = p, 
device = 'pdf', width = 14, height = 10, useDingbats = FALSE)


#And Tac1 and cholinergic coexpression
genes <- c('Tac1',"Chat","Slc5a7")
expr_mat <- assay(all_donor_sce, "logcounts")

cluster_to_annotate = "meta_clust_celltype_annot"
clusters <- unique(colData(all_donor_sce)[[cluster_to_annotate]])

coexp_df <- lapply(clusters, function(cl) {
  cells <- colData(all_donor_sce)[[cluster_to_annotate]] == cl
  mat_sub <- expr_mat[, cells, drop = FALSE]
  pairwise_coexpression(mat_sub, genes, cluster_name = cl)
}) %>%
  bind_rows()

coexp_df$gene1 <- factor(coexp_df$gene1, levels = genes)
coexp_df$gene2 <- factor(coexp_df$gene2, levels = rev(genes))


#Edited the original code to have grey be between 0-5%, previously the very low percentages were difficult to see.
p <- ggplot(coexp_df, aes(x = gene1, y = gene2, fill = percent)) +
  geom_tile(color = "grey70", linewidth = 0.3) +
  facet_wrap(~ cluster, nrow = 5) +
  scale_fill_gradientn(
    colours = c("grey95","#f1e2c6", "#f1e2c6", "#f0c94a", "#df8b27", "#d92523", "#8b0d19"),
    values = c(0, 0.05, 0.20, 0.40, 0.60, 0.8, 1),
    limits = c(0, 100),
    breaks = c( 5, 20, 40, 60, 80, 100),
    name = "Percent of cells expressing\ntwo genes"
  ) +
  coord_equal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  xlab(NULL) +
  ylab(NULL) + ggtitle("Co-expression of SubP and Chol markers in Wallace mouse dataset")

print(p)

ggsave(path = plot_path, filename = 'SubP_Chol_coexpression_Wallace_mouse.pdf', plot = p, 
device = 'pdf', width = 14, height = 10, useDingbats = FALSE)







