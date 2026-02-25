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


#I think it will be clearer to just focus on the habenula data, so filter the non-neurons in Wallace and use the neuronal subset from Hashikawa
wallace_sce = wallace_sce[, !wallace_sce$celltype %in% c('Astrocytes','Polydendrocytes','Differentiating Oligodendrocytes','Oligodendrocytes',
'Pericytes','Fibroblasts','Endothelial','Macrophages', 'Microglia')]
table(wallace_sce$celltype)

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


studies <- unique(hashikawa_sce_sub$stim)
hashikawa_sce_sub_list <- lapply(studies, function(study) {
  hashikawa_sce_sub[, hashikawa_sce_sub$stim == study]
})
names(hashikawa_sce_sub_list) <- studies


########################
#MetaNeighbor
#########################

#Get single SCE object
all_donor_sce = mergeSCE(c(wallace_sce_list, hashikawa_sce_sub_list))
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


#From these matches
#Hashikawa
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













