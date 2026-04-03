##################################
#Compare Wallace to the Hashikawa dataset
#path is /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data
#First check out what is there, probably best to just use what has been processed and annotated already
###################################

library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
library(MetaNeighbor)
library(dplyr)
library(ggplot2)
library(here)

here::here()


#Path to save any generated data
new_data_path = here('processed-data', '98_external_Hb_comparisons', '04_wallace_hashikawa_mouse')
#Path to already generated data
wallace_data_path = here('processed-data', '98_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')
wallace_meta_path = here('processed-data', '98_external_Hb_comparisons', '03_metaMarkers_wallace_2019')
#Path to plot directory
plot_path = here('plots', '98_external_Hb_comparisons', '04_wallace_hashikawa_mouse')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


#Source the bubble plot functions
source(here('code','98_external_Hb_comparisons', 'bubble_plot_functions.R'))



#From the /dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/code/09_cross_species_analysis/03_Hashikawa_data_prep.R
#What I want is the sce_mouse_habenula.Rdata object

hashikawa_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data'
list.files(hashikawa_path)

#Loads as sce_mouse_sub, list of sce objects
load(paste0(hashikawa_path, '/sce_mouse_habenula.Rdata'))
hashikawa_sce = sce_mouse_sub$all
hashikawa_sce_sub = sce_mouse_sub$neuron
rm(sce_mouse_sub)

#Swap the rownames to gene symbols
rowData(hashikawa_sce)
rownames(hashikawa_sce) = rowData(hashikawa_sce)$Symbol
rownames(hashikawa_sce_sub) = rowData(hashikawa_sce_sub)$Symbol

#What do the annotations look like
#Only has Neuron1, 2, 3, etc labels
#Is this only lateral habenula?
colnames(colData(hashikawa_sce))
table(hashikawa_sce$celltype)
table(hashikawa_sce_sub$celltype)
#Probably split by control and stimulus
table(hashikawa_sce$stim)
table(hashikawa_sce_sub$stim)

#Load up the Wallace data
wallace_sce = readRDS(paste0(
  wallace_data_path, '/all_donor_sce_with_denovo_clusters.rds'))

#add in the annotated metdata
current_mouse_metadata = readRDS(paste0(wallace_meta_path, '/wallace_mouse_metaclust_celltype_annot_metadata.rds'))
colData(wallace_sce) = S4Vectors::DataFrame(current_mouse_metadata)

wallace_sce


#Add cpm to the hashikawa data
assay(hashikawa_sce, "cpm") = MetaMarkers::convert_to_cpm(assay(hashikawa_sce, "counts"))
assay(hashikawa_sce_sub, "cpm") = MetaMarkers::convert_to_cpm(assay(hashikawa_sce_sub, "counts"))

#Filter for the shared present genes
present_genes = intersect(intersect(rownames(wallace_sce), rownames(hashikawa_sce)), rownames(hashikawa_sce_sub))
length(present_genes)

hashikawa_sce = hashikawa_sce[present_genes, ]
hashikawa_sce_sub = hashikawa_sce_sub[present_genes, ]
wallace_sce = wallace_sce[present_genes, ]

#Sanity check
table(rownames(hashikawa_sce) %in% rownames(wallace_sce))
table(rownames(hashikawa_sce_sub) %in% rownames(wallace_sce))

#Rename cell-type metadata to match
wallace_sce$celltype = wallace_sce$meta_clust_celltype_annot


#Do a comparison including the non-neurons, and then just the neurons
wallace_neuron_sce = wallace_sce[, !wallace_sce$celltype %in% c('Astrocytes','Polydendrocytes','Differentiating Oligodendrocytes','Oligodendrocytes',
'Pericytes','Fibroblasts','Endothelial','Macrophages', 'Microglia')]

table(wallace_neuron_sce$celltype)
table(wallace_sce$celltype)
table(hashikawa_sce$celltype)



########
#Do the comparison including the non-neurons first
########

#Split wallace by donor and Hashikawa by simulus condition
studies <- unique(wallace_sce$study_id)
wallace_sce_list <- lapply(studies, function(study) {
  wallace_sce[, wallace_sce$study_id == study]
})
names(wallace_sce_list) <- studies


studies <- unique(hashikawa_sce$stim)
hashikawa_sce_list <- lapply(studies, function(study) {
  hashikawa_sce[, hashikawa_sce$stim == study]
})
names(hashikawa_sce_list) <- studies



########################
#MetaNeighbor
#########################

#Get single SCE object
all_donor_sce = mergeSCE(c(wallace_sce_list, hashikawa_sce_list))
View(as.data.frame(colData(all_donor_sce)))

#Ignore the outlier cells
all_donor_sce = all_donor_sce[,
  all_donor_sce$celltype != 'outliers'
]


#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(
  dat = all_donor_sce,
  min_recurrence = 2,
  exp_labels = all_donor_sce$study_id
)
length(global_hvgs)
keep_global_hvgs = global_hvgs[1:2000]


MN_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$celltype,
  fast_version = TRUE
)

#Plot allby-all AUROC heatmap
plotHeatmap(
  MN_aurocs,
  show_dendro = TRUE,
  show_labels = TRUE,
  cex = .5
)
title("MetaNeighbor Mouse Wallace vs Hashikawa: 2000 HVGs")

pdf(paste0(plot_path, '/all_vs_all_MN_aurocs_heatmap.pdf'), width = 10, height = 8)
plotHeatmap(
  MN_aurocs,
  show_dendro = TRUE,
  show_labels = TRUE,
  cex = .5
)
title("MetaNeighbor Mouse Wallace vs Hashikawa: 2000 HVGs")
dev.off()


#And the best versus next approach

MN_best_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$celltype,
  fast_version = TRUE,
  one_vs_best = TRUE,
  symmetric_output = FALSE
)

#Plot best_vs_next AUROC heatmap
plotHeatmap(
  MN_best_aurocs,
  show_dendro = TRUE,
  show_labels = TRUE,
  cex = .5
)
title("MetaNeighbor BvsNext Mouse Wallace vs Hashikawa: 2000 HVGs")

pdf(paste0(plot_path, '/best_vs_next_MN_aurocs_heatmap.pdf'), width = 10, height = 8)
plotHeatmap(
  MN_best_aurocs,
  show_dendro = TRUE,
  show_labels = TRUE,
  cex = .5
)
title("MetaNeighbor BvsNext Mouse Wallace vs Hashikawa: 2000 HVGs")
dev.off()


cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$study_id, all_donor_sce$celltype, size_factor = 3)

pdf(paste0(plot_path, '/MN_cluster_graph.pdf'), width = 10, height = 8)
cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$study_id, all_donor_sce$celltype, size_factor = 3)
dev.off()


table(hashikawa_sce$celltype)

#From these matches
#Hashikawa neurons
#Neuron1 - Cholinergic
#Neuron2 - LHb1
#Neuron3 - SubP-Cholinergic
#Neuron4 - LHb1
#Neuron5 - LHb2
#Neuron6 - Cholinergic
#Neuron7 - Cholinergic
#Neuron8 - LHb1

rownames(hashikawa_sce) = toupper(rownames(hashikawa_sce))

p_bubble = get_bubble_plot_sce(hashikawa_sce, 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1', 'TAC2','GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1'), sample_name = "Hashikawa mouse Habenula", group_col = "celltype")

p_bubble[[1]]
p_bubble[[2]]

pdf(paste0(plot_path, '/Hashikawa_Hab_marker_bubbles_meanExp.pdf'), width = 10, height = 8)
p_bubble[[1]]
dev.off()
pdf(paste0(plot_path, '/Hashikawa_Hab_marker_bubbles_Zscore_meanExp.pdf'), width = 10, height = 8)
p_bubble[[2]]
dev.off()


#Altogether, very clear matches between the mouse datasets.
#Substance P and cholinergic populations clear

rownames(wallace_sce) = toupper(rownames(wallace_sce))
p_bubble = get_bubble_plot_sce(wallace_sce[ , wallace_sce$celltype != 'outliers'], 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1', 'TAC2','GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1'), sample_name = "Wallace mouse Habenula", group_col = "celltype")

p_bubble[[1]]
p_bubble[[2]]


pdf(paste0(plot_path, '/Wallace_Hab_marker_bubbles_meanExp.pdf'), width = 10, height = 8)
p_bubble[[1]]
dev.off()
pdf(paste0(plot_path, '/Wallace_Hab_marker_bubbles_Zscore_meanExp.pdf'), width = 10, height = 8)
p_bubble[[2]]
dev.off()

#The inhibitory markers are still expressed in Wallace LHb_1 over LHb_2, seems consistent across the datasets
#The lack of VGLUT1 expression in all lateral Hab clusters is really clear across both datasets


#Save metadata for the Hashikawa dataset with annotations matched to the wallace dataset
table(hashikawa_sce$celltype)
#metacluster annotations from all the above
meta_annot_vec = c('MHb_cholinergic.1' = 'Neuron1',
                   'LHb_1.1' = 'Neuron2',
                   'MHb_subP_cholinergic' = 'Neuron3',
                   'LHb_1.2' = 'Neuron4',
                   'LHb_2' = 'Neuron5',
                   'MHb_cholinergic.2' = 'Neuron6',
                   'MHb_cholinergic.3' = 'Neuron7',
                   'LHb_1.3' = 'Neuron8',
                   'Astrocyte1' = 'Astrocyte1',
                   'Astrocyte2' = 'Astrocyte2',
                   'Endothelial' = 'Endothelial',
                   'Epen' = 'Epen',
                   'Microglia' = 'Microglia',
                   'Mural' = 'Mural',
                   'Oligo1' = 'Oligo1',
                   'Oligo2' = 'Oligo2',
                   'Oligo3' = 'Oligo3',
                   'OPC1' = 'OPC1',
                   'OPC2' = 'OPC2',
                   'OPC3' = 'OPC3'
)


meta_annot_vec  = setNames(names(meta_annot_vec), meta_annot_vec)


hashikawa_sce$meta_clust_celltype_annot = unname(meta_annot_vec[hashikawa_sce$celltype])
table(hashikawa_sce$meta_clust_celltype_annot)

#Save the metadata as a data.frame to add to the seurat data object later
full_hashikawa_metadata = colData(hashikawa_sce)
saveRDS(full_hashikawa_metadata, paste0(new_data_path, '/hashikawa_mouse_metaclust_celltype_annot_metadata.rds'))






#And now with just the neuronal subset

#Split wallace by donor and Hashikawa by simulus condition
studies <- unique(wallace_neuron_sce$study_id)
wallace_neuron_sce_list <- lapply(studies, function(study) {
  wallace_neuron_sce[, wallace_neuron_sce$study_id == study]
})
names(wallace_neuron_sce_list) <- studies


studies <- unique(hashikawa_sce_sub$stim)
hashikawa_sce_sub_list <- lapply(studies, function(study) {
  hashikawa_sce_sub[, hashikawa_sce_sub$stim == study]
})
names(hashikawa_sce_sub_list) <- studies



########################
#MetaNeighbor
#########################

#Get single SCE object
all_donor_sce = mergeSCE(c(wallace_neuron_sce_list, hashikawa_sce_sub_list))
View(as.data.frame(colData(all_donor_sce)))

#Ignore the outlier cells
all_donor_sce = all_donor_sce[,
  all_donor_sce$celltype != 'outliers'
]


#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(
  dat = all_donor_sce,
  min_recurrence = 2,
  exp_labels = all_donor_sce$study_id
)
length(global_hvgs)
keep_global_hvgs = global_hvgs[1:2000]


MN_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$celltype,
  fast_version = TRUE
)

#Plot allby-all AUROC heatmap
plotHeatmap(
  MN_aurocs,
  show_dendro = TRUE,
  show_labels = TRUE,
  cex = .5
)
title("MetaNeighbor Mouse Wallace vs Hashikawa: 2000 HVGs")

pdf(paste0(plot_path, '/all_vs_all_MN_just_neurons_aurocs_heatmap.pdf'), width = 10, height = 8)
plotHeatmap(
  MN_aurocs,
  show_dendro = TRUE,
  show_labels = TRUE,
  cex = .5
)
title("MetaNeighbor Mouse Wallace vs Hashikawa: 2000 HVGs")
dev.off()


#And the best versus next approach

MN_best_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$celltype,
  fast_version = TRUE,
  one_vs_best = TRUE,
  symmetric_output = FALSE
)

#Plot best_vs_next AUROC heatmap
plotHeatmap(
  MN_best_aurocs,
  show_dendro = TRUE,
  show_labels = TRUE,
  cex = .5
)
title("MetaNeighbor BvsNext Mouse Wallace vs Hashikawa: 2000 HVGs")

pdf(paste0(plot_path, '/best_vs_next_MN_just_neurons_aurocs_heatmap.pdf'), width = 10, height = 8)
plotHeatmap(
  MN_best_aurocs,
  show_dendro = TRUE,
  show_labels = TRUE,
  cex = .5
)
title("MetaNeighbor BvsNext Mouse Wallace vs Hashikawa: 2000 HVGs")
dev.off()


cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$study_id, all_donor_sce$celltype, size_factor = 3)

pdf(paste0(plot_path, '/MN_cluster_graph_just_neurons.pdf'), width = 10, height = 8)
cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$study_id, all_donor_sce$celltype, size_factor = 3)
dev.off()


table(hashikawa_sce$celltype)




#MHb2, 3, 4, 6 : Cholinergic
#MHb 5 : Substance P
#MHb1: subP_cholinergic

#LHb1, 2, 3, 4, 6: Wallace LHb1
#LHb5: Wallace LHb2

#Check out the markers in the hashikawa dataset


rownames(hashikawa_sce_sub) = toupper(rownames(hashikawa_sce_sub))

p_bubble = get_bubble_plot_sce(hashikawa_sce_sub, 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1', 'TAC2','GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1'), sample_name = "Hashikawa mouse Habenula", group_col = "celltype")

p_bubble[[1]]
p_bubble[[2]]

pdf(paste0(plot_path, '/Hashikawa_Hab_just_neurons_marker_bubbles_meanExp.pdf'), width = 10, height = 8)
p_bubble[[1]]
dev.off()
pdf(paste0(plot_path, '/Hashikawa_Hab_just_neurons_marker_bubbles_Zscore_meanExp.pdf'), width = 10, height = 8)
p_bubble[[2]]
dev.off()


#Altogether, very clear matches between the mouse datasets.
#Substance P and cholinergic populations clear

rownames(wallace_neuron_sce) = toupper(rownames(wallace_neuron_sce))
p_bubble = get_bubble_plot_sce(wallace_neuron_sce[ , wallace_neuron_sce$celltype != 'outliers'], 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1', 'TAC2','GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1'), sample_name = "Wallace mouse Habenula", group_col = "celltype")

p_bubble[[1]]
p_bubble[[2]]


pdf(paste0(plot_path, '/Wallace_Hab_just_neurons_marker_bubbles_meanExp.pdf'), width = 10, height = 8)
p_bubble[[1]]
dev.off()
pdf(paste0(plot_path, '/Wallace_Hab_just_neurons_marker_bubbles_Zscore_meanExp.pdf'), width = 10, height = 8)
p_bubble[[2]]
dev.off()

#The inhibitory markers are still expressed in Wallace LHb_1 over LHb_2, seems consistent across the datasets
#The lack of VGLUT1 expression in all lateral Hab clusters is really clear across both datasets


#Save metadata for the Hashikawa dataset with annotations matched to the wallace dataset

#metacluster annotations from all the above
meta_annot_vec = c('MHb_subP_cholinergic' = 'MHb1',
                   'MHb_cholinergic_1' = 'MHb2',
                   'MHb_cholinergic_2' = 'MHb3',
                   'MHb_cholinergic_3' = 'MHb4',
                   'MHb_subP' = 'MHb5',
                   'MHb_cholinergic_4' = 'MHb6',
                   'LHb_1_1' = 'LHb1',
                   'LHb_1_2' = 'LHb2',
                   'LHb_1_3' = 'LHb3',
                   'LHb_1_4' = 'LHb4',
                   'LHb_2' = 'LHb5',
                   'LHb_1_5' = 'LHb6'
)

meta_annot_vec  = setNames(names(meta_annot_vec), meta_annot_vec)


hashikawa_sce_sub$meta_clust_celltype_annot = unname(meta_annot_vec[hashikawa_sce_sub$celltype])
table(hashikawa_sce_sub$meta_clust_celltype_annot)

#Save the metadata as a data.frame to add to the seurat data object later
full_hashikawa_metadata = colData(hashikawa_sce_sub)
saveRDS(full_hashikawa_metadata, paste0(new_data_path, '/hashikawa_mouse_neuron_metaclust_celltype_annot_metadata.rds'))

##########################
#Look for potential lateral habenula GABAergic cells
#########################

#Switch the full dataset to Seurat for the bubble plots
hashikawa_seurat = as.Seurat(hashikawa_sce_sub, counts = "counts", data = "cpm")

#HVGs
hashikawa_seurat <- FindVariableFeatures(hashikawa_seurat, selection.method = "vst", nfeatures = 2000)

#Scale data
all.genes <- rownames(hashikawa_seurat)
hashikawa_seurat <- ScaleData(hashikawa_seurat, features = all.genes)

#PCA
hashikawa_seurat  <- RunPCA(hashikawa_seurat , features = VariableFeatures(object = hashikawa_seurat ))
DimPlot(hashikawa_seurat, reduction = "pca") + NoLegend()

#UMAP
hashikawa_seurat  <- RunUMAP(hashikawa_seurat , dims = 1:20)
DimPlot(hashikawa_seurat, reduction = "umap", group.by = 'meta_clust_celltype_annot')
DimPlot(hashikawa_seurat, reduction = "umap", group.by = 'stim')


custom_markers = c('Gad1', 'Gad2', 'Slc32a1','Pvalb', 'Slc17a6', 'Slc17a7')

p_bubble_gaba_glut_markers = get_bubble_plot(hashikawa_seurat , custom_markers, 'Hashikawa: all mouse samples', 
group_col = 'meta_clust_celltype_annot')
p_bubble_gaba_glut_markers


#The LHb1 has the highest Gad2 expression and VGAT expression, but both are low in general and small percentages. This matches the trend in humans
mouse_gad1_p = FeaturePlot(hashikawa_seurat, features = "Gad1", reduction = "umap", pt.size = 1, slot = 'scale.data') +
  scale_color_gradient(low = "white", high = "red", name = 'z-score')

mouse_gad2_p = FeaturePlot(hashikawa_seurat, features = "Gad2", reduction = "umap", pt.size = 1, slot = 'scale.data') +
  scale_color_gradient(low = "white", high = "red", name = 'z-score')

mouse_pvalb_p = FeaturePlot(hashikawa_seurat, features = "Pvalb", reduction = "umap", pt.size = 1, slot = 'scale.data') +
  scale_color_gradient(low = "white", high = "red", name = 'z-score')

mouse_vgat_p = FeaturePlot(hashikawa_seurat, features = "Slc32a1", reduction = "umap", pt.size = 1, slot = 'scale.data') +
  scale_color_gradient(low = "white", high = "red", name = 'z-score')


#Look at the co-expression of specific genes
# Get expression data
umap_data <- as.data.frame(Embeddings(hashikawa_seurat, reduction = "umap"))
umap_data$GAD2 <- FetchData(hashikawa_seurat, vars = "Gad2", slot = "data")[, 1]
umap_data$PVALB <- FetchData(hashikawa_seurat, vars = "Pvalb", slot = "data")[, 1]
umap_data$SLC17A6 <- FetchData(hashikawa_seurat, vars = "Slc17a6", slot = "data")[, 1]

# Create a coexpression category
umap_data$coexpression <- ifelse(umap_data$GAD2 > 0 & umap_data$SLC17A6 > 0, "Both",
                                  ifelse(umap_data$GAD2 > 0, "Gad2 only",
                                         ifelse(umap_data$SLC17A6 > 0, "Slc17a6 only", "Neither")))

gad2_vglut2_p = ggplot(umap_data, aes(x = umap_1, y = umap_2, color = coexpression)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Both" = "purple", "Gad2 only" = "red", "Slc17a6 only" = "blue", "Neither" = "lightgrey")) +
  theme_bw() +
  labs(title = "Gad2 and Slc17a6 Co-expression: Hashikawa mouse")


umap_data$coexpression <- ifelse(umap_data$GAD2 > 0 & umap_data$PVALB > 0, "Both",
                                  ifelse(umap_data$GAD2 > 0, "Gad2 only",
                                         ifelse(umap_data$PVALB > 0, "Pvalb only", "Neither")))

gad2_pvalb_p = ggplot(umap_data, aes(x = umap_1, y = umap_2, color = coexpression)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Both" = "purple", "Gad2 only" = "red", "Pvalb only" = "blue", "Neither" = "lightgrey")) +
  theme_bw() +
  labs(title = "Gad2 and Pvalb Co-expression: Hashikawa mouse")



p_bubble_gaba_glut_markers[[1]]
p_bubble_gaba_glut_markers[[2]]

mouse_gad1_p 
mouse_gad2_p
mouse_pvalb_p
mouse_vgat_p

gad2_vglut2_p 
gad2_pvalb_p 

ggsave(p_bubble_gaba_glut_markers[[1]], filename = 'hashikawa_mouse_gaba_glut_meta_annots_bubble.pdf', path = plot_path,
device = 'pdf', width = 10, height = 8)

ggsave(p_bubble_gaba_glut_markers[[2]], filename = 'hashikawa_mouse_gaba_glut_zscore_meta_annots_bubble.pdf', path = plot_path,
device = 'pdf', width = 10, height = 8)

ggsave(mouse_gad1_p , filename = 'hashikawa_mouse_gad1_exp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(mouse_pvalb_p , filename = 'hashikawa_mouse_pvalb_exp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(mouse_gad2_p , filename = 'hashikawa_mouse_gad2_exp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(mouse_vgat_p , filename = 'hashikawa_mouse_vgat_exp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(gad2_vglut2_p , filename = 'hashikawa_mouse_gad2_vglut2_coexp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(gad2_pvalb_p , filename = 'hashikawa_mouse_gad2_pvalb_coexp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

#Coexpression plots of the GABA-Glut genes


logcounts(hashikawa_sce_sub) = log1p(assay(hashikawa_sce_sub, "cpm"))
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
expr_mat <- assay(hashikawa_sce_sub, "logcounts")

cluster_to_annotate = "meta_clust_celltype_annot"
clusters <- unique(colData(hashikawa_sce_sub)[[cluster_to_annotate]])

coexp_df <- lapply(clusters, function(cl) {
  cells <- colData(hashikawa_sce_sub)[[cluster_to_annotate]] == cl
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
  ylab(NULL) + ggtitle("Co-expression of GABA and Glut markers in Hashikawa mouse dataset")

print(p)

ggsave(path = plot_path, filename = 'GABA_Glut_coexpression_Hashikawa_mouse.pdf', plot = p, 
device = 'pdf', width = 14, height = 10, useDingbats = FALSE)

#And Tac1 and cholinergic coexpression
genes <- c('Tac1',"Chat","Slc5a7")
expr_mat <- assay(hashikawa_sce_sub, "logcounts")

cluster_to_annotate = "meta_clust_celltype_annot"
clusters <- unique(colData(hashikawa_sce_sub)[[cluster_to_annotate]])

coexp_df <- lapply(clusters, function(cl) {
  cells <- colData(hashikawa_sce_sub)[[cluster_to_annotate]] == cl
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
  ylab(NULL) + ggtitle("Co-expression of SubP and Chol markers in Hashikawa mouse dataset")

print(p)

ggsave(path = plot_path, filename = 'SubP_Chol_coexpression_Hashikawa_mouse.pdf', plot = p, 
device = 'pdf', width = 14, height = 10, useDingbats = FALSE)


#Compute mouse meta-markers at a broad annotation level
#For this purpose, just group all the non-neurons together
#Group the hashikawa data to match the wallace annotations
table(hashikawa_sce$meta_clust_celltype_annot)

broad_annot_vec = c('MHb_cholinergic' = 'MHb_cholinergic.1',
                   'MHb_cholinergic' = 'MHb_cholinergic.2',
                   'MHb_cholinergic' = 'MHb_cholinergic.3',
                   'MHb_subP_cholinergic' = 'MHb_subP_cholinergic',
                   'LHb_1' = 'LHb_1.1',
                   'LHb_1' = 'LHb_1.2',
                   'LHb_1' = 'LHb_1.3',
                   'LHb_2' = 'LHb_2',
                   'non_neurons' = 'Astrocyte1',
                   'non_neurons' = 'Astrocyte2',
                   'non_neurons' = 'Endothelial',
                   'non_neurons' = 'Epen',
                   'non_neurons' = 'Microglia',
                   'non_neurons' = 'Mural',
                   'non_neurons' = 'Oligo1',
                   'non_neurons' = 'Oligo2',
                   'non_neurons' = 'Oligo3',
                   'non_neurons' = 'OPC1',
                   'non_neurons' = 'OPC2',
                   'non_neurons' = 'OPC3'
)

broad_annot_vec  = setNames(names(broad_annot_vec), broad_annot_vec)
hashikawa_sce$broad_celltype_annot = unname(broad_annot_vec[hashikawa_sce$meta_clust_celltype_annot])
table(hashikawa_sce$broad_celltype_annot, hashikawa_sce$meta_clust_celltype_annot )


table(wallace_sce$meta_clust_celltype_annot)
broad_annot_vec = c('MHb_cholinergic' = 'MHb_cholinergic',
                   'MHb_subP' = 'MHb_subP',
                   'MHb_subP_cholinergic' = 'MHb_subP_cholinergic',
                   'LHb_1' = 'LHb_1',
                   'LHb_2' = 'LHb_2',
                   'non_neurons' = 'Astrocytes',
                   'non_neurons' = 'Differentiating Oligodendrocytes',
                   'non_neurons' = 'Endothelial',
                   'non_neurons' = 'Fibroblasts',
                   'non_neurons' = 'Macrophages',
                   'non_neurons' = 'Microglia',
                   'non_neurons' = 'Oligodendrocytes',
                   'non_neurons' = 'Pericytes',
                   'non_neurons' = 'Polydendrocytes',
                   'outliers' = 'outliers'
)

broad_annot_vec  = setNames(names(broad_annot_vec), broad_annot_vec)
wallace_sce$broad_celltype_annot = unname(broad_annot_vec[wallace_sce$meta_clust_celltype_annot])
table(wallace_sce$broad_celltype_annot, wallace_sce$meta_clust_celltype_annot )


#And the neuronal subsets

table(hashikawa_sce_sub$meta_clust_celltype_annot)

broad_annot_vec = c('MHb_cholinergic' = 'MHb_cholinergic_1',
                   'MHb_cholinergic' = 'MHb_cholinergic_2',
                   'MHb_cholinergic' = 'MHb_cholinergic_3',
                   'MHb_cholinergic' = 'MHb_cholinergic_4',
                   'MHb_subP_cholinergic' = 'MHb_subP_cholinergic',
                   'MHb_subP' = 'MHb_subP',
                   'LHb_1' = 'LHb_1_1',
                   'LHb_1' = 'LHb_1_2',
                   'LHb_1' = 'LHb_1_3',
                   'LHb_1' = 'LHb_1_4',
                   'LHb_1' = 'LHb_1_5',
                   'LHb_2' = 'LHb_2'
)

broad_annot_vec  = setNames(names(broad_annot_vec), broad_annot_vec)
hashikawa_sce_sub$broad_celltype_annot = unname(broad_annot_vec[hashikawa_sce_sub$meta_clust_celltype_annot])
table(hashikawa_sce_sub$broad_celltype_annot, hashikawa_sce_sub$meta_clust_celltype_annot )

#Dont need to group for wallace
table(wallace_neuron_sce$meta_clust_celltype_annot)
wallace_neuron_sce$broad_celltype_annot = wallace_neuron_sce$meta_clust_celltype_annot




#Get markers and mouse meta-markers
hashikawa_cntl_all = hashikawa_sce[ ,hashikawa_sce$stim == 'cntl' ]
hashikawa_stim_all = hashikawa_sce[ ,hashikawa_sce$stim == 'stim' ]

wallace_d1_all = wallace_sce[ , wallace_sce$study_id == 'hab_160822' ]
wallace_d2_all = wallace_sce[ , wallace_sce$study_id == 'hab_161102' ]
wallace_d3_all = wallace_sce[ , wallace_sce$study_id == 'hab_161103' ]
wallace_d4_all = wallace_sce[ , wallace_sce$study_id == 'hab_161105' ]

all_markers = list(
  hashikawa_cntl_all = compute_markers(assay(hashikawa_cntl_all, "cpm"), hashikawa_cntl_all$broad_celltype_annot),
  hashikawa_stim_all = compute_markers(assay(hashikawa_stim_all, "cpm"), hashikawa_stim_all$broad_celltype_annot),
  wallace_d1_all = compute_markers(assay(wallace_d1_all, "cpm"), wallace_d1_all$broad_celltype_annot),
  wallace_d2_all = compute_markers(assay(wallace_d2_all, "cpm"), wallace_d2_all$broad_celltype_annot),
  wallace_d3_all = compute_markers(assay(wallace_d3_all, "cpm"), wallace_d3_all$broad_celltype_annot),
  wallace_d4_all = compute_markers(assay(wallace_d4_all, "cpm"), wallace_d4_all$broad_celltype_annot)
)

mouse_all_meta_markers = make_meta_markers(all_markers, detailed_stats = TRUE)

mouse_all_meta_markers %>% filter(rank <= 10) %>% View()

export_meta_markers(mouse_all_meta_markers, paste0(new_data_path, '/mouse_hab_all_celltypes_meta_markers.csv'), names(all_markers))



#Just the neuron subsets
hashikawa_cntl_all = hashikawa_sce_sub[ ,hashikawa_sce_sub$stim == 'cntl' ]
hashikawa_stim_all = hashikawa_sce_sub[ ,hashikawa_sce_sub$stim == 'stim' ]

wallace_d1_all = wallace_neuron_sce[ , wallace_neuron_sce$study_id == 'hab_160822' ]
wallace_d2_all = wallace_neuron_sce[ , wallace_neuron_sce$study_id == 'hab_161102' ]
wallace_d3_all = wallace_neuron_sce[ , wallace_neuron_sce$study_id == 'hab_161103' ]
wallace_d4_all = wallace_neuron_sce[ , wallace_neuron_sce$study_id == 'hab_161105' ]

neuron_markers = list(
  hashikawa_cntl_all = compute_markers(assay(hashikawa_cntl_all, "cpm"), hashikawa_cntl_all$broad_celltype_annot),
  hashikawa_stim_all = compute_markers(assay(hashikawa_stim_all, "cpm"), hashikawa_stim_all$broad_celltype_annot),
  wallace_d1_all = compute_markers(assay(wallace_d1_all, "cpm"), wallace_d1_all$broad_celltype_annot),
  wallace_d2_all = compute_markers(assay(wallace_d2_all, "cpm"), wallace_d2_all$broad_celltype_annot),
  wallace_d3_all = compute_markers(assay(wallace_d3_all, "cpm"), wallace_d3_all$broad_celltype_annot),
  wallace_d4_all = compute_markers(assay(wallace_d4_all, "cpm"), wallace_d4_all$broad_celltype_annot)
)

mouse_neuron_meta_markers = make_meta_markers(neuron_markers, detailed_stats = TRUE)

mouse_neuron_meta_markers %>% filter(rank <= 25) %>% View()

export_meta_markers(mouse_neuron_meta_markers, paste0(new_data_path, '/mouse_hab_neuron_meta_markers.csv'), names(all_markers))


mouse_neuron_meta_markers = read_meta_markers(paste0(new_data_path, '/mouse_hab_neuron_meta_markers.csv.gz'))
mouse_neuron_meta_markers %>% filter(rank <= 25) %>% View()


mouse_all_meta_markers = read_meta_markers(paste0(new_data_path, '/mouse_hab_all_celltypes_meta_markers.csv.gz'))
mouse_all_meta_markers %>% filter(rank <= 25) %>% View()
