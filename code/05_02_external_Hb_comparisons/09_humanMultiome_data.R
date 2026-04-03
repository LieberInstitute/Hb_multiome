#Incorportating the human multiome habenula data from the atlasing project. Compare the current transcriptional clusters to the zebrafish, mouse, and Yalcinbas human data

library(SingleCellExperiment)
library(Seurat)
library(MetaNeighbor)
library(MetaMarkers)
library(qs2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(here)

here::here()


#Path to the human multiome data
#Using the v5.rds object, says seurat, but is SCE object
multiome_path = here('processed-data', '05_01_drop_doublets','01_drop_doublets_and_reDimReduce')
list.files(multiome_path)

#Path to the Yalcinbas pilot data
#Going with the official_final_sce.RDATA
yalcinbas_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/sce_objects'
list.files(yalcinbas_path)

#Nonhuman data paths
zeb_data_path = here('processed-data', '05_02_external_Hb_comparisons','07_cross_species_hab')
mouse_data_path = here('processed-data', '05_02_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')
wallace_03_data_path = here('processed-data', '05_02_external_Hb_comparisons', '03_metaMarkers_wallace_2019')

hashikawa_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data'
hashikawa_data_path = here('processed-data','05_02_external_Hb_comparisons','04_wallace_hashikawa_mouse')

#Path to orthologs
path_to_orthologs = here('processed-data', '05_02_external_Hb_comparisons', 'human_mouse_zebrafish_orthologs.txt.gz')


#Path for new data generated
new_data_path = here('processed-data', '05_02_external_Hb_comparisons','09_humanMultiome_data')
#Path to plot directory
plot_path = here('plots', '05_02_external_Hb_comparisons', '09_humanMultiome_data')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)

#Source the bubble plot functions
source(here('code','05_02_external_Hb_comparisons', 'bubble_plot_functions.R'))


#Multiome human data
multiome_sce = qs_read(paste0(multiome_path, '/reprocessed_doubletRemoved_multiomeHab_SCE.qs2'))
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))

names(colData(multiome_sce))
#Low resolution clusters
table(multiome_sce$merged_cluster)
#Resolution used in other comparisons
table(multiome_sce$mid_cluster)
#Donor ID
table(multiome_sce$orig.ident)

#Yalcinbas data
#Loads as an object labeled 'sce'
#16437 cells
load(paste0(yalcinbas_path, '/official_final_sce.RDATA'))
yalcinbas_sce = sce
rm(sce)

table(yalcinbas_sce$final_Annotations)

assay(yalcinbas_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(yalcinbas_sce, 'counts'))



#Load up the zebrafish and mouse data, filter down to the present genes
#Zebrafish data, using the summed paralog version
zeb_sce = readRDS(file = paste0(zeb_data_path, '/adult_zebrafish_summed_paralogs.rds'))

#Wallace mouse data
all_mouse_sce = readRDS(paste0(mouse_data_path, '/all_donor_sce_with_denovo_clusters.rds'))

#Hashikawa mouse data
load(paste0(hashikawa_path, '/sce_mouse_habenula.Rdata'))
hashikawa_sce_sub = sce_mouse_sub$neuron #Using just the neuron subset
rm(sce_mouse_sub)
#Gene symbols as the rownames
rownames(hashikawa_sce_sub) = rowData(hashikawa_sce_sub)$Symbol
#Add CPM
assay(hashikawa_sce_sub, "cpm") = MetaMarkers::convert_to_cpm(assay(hashikawa_sce_sub, "counts"))


hu_mu_zf_ortholog_df = data.table::fread(path_to_orthologs)
table(hu_mu_zf_ortholog_df$`Mouse homology type`)

#Get the 1to1 orthologs for the mouse
hu_mu_zf_ortholog_df <- hu_mu_zf_ortholog_df %>%
  filter(
    `Mouse homology type` == 'ortholog_one2one'
  ) %>%
  filter(!duplicated(`Gene name`))
dim(hu_mu_zf_ortholog_df)

index = match(rownames(all_mouse_sce), hu_mu_zf_ortholog_df$`Mouse gene name`)
rownames(all_mouse_sce) = hu_mu_zf_ortholog_df$`Gene name`[index]

index = match(rownames(hashikawa_sce_sub), hu_mu_zf_ortholog_df$`Mouse gene name`)
rownames(hashikawa_sce_sub) = hu_mu_zf_ortholog_df$`Gene name`[index]

#Get the set of genes that are present in all datasets, and are 1-to-1 orthologs
#9522 across all three, which is pretty good
all_present_genes = Reduce(intersect, 
  list(rownames(all_mouse_sce), rownames(hashikawa_sce_sub), rownames(zeb_sce), rownames(yalcinbas_sce), rownames(multiome_sce)))
all_present_genes = all_present_genes[!is.na(all_present_genes)]
length(all_present_genes)


#Filter the SCE objects for the present genes
zeb_sce = zeb_sce[rownames(zeb_sce) %in% all_present_genes, ]
all_mouse_sce = all_mouse_sce[rownames(all_mouse_sce) %in% all_present_genes, ]
hashikawa_sce_sub = hashikawa_sce_sub[rownames(hashikawa_sce_sub) %in% all_present_genes, ]
yalcinbas_sce = yalcinbas_sce[rownames(yalcinbas_sce) %in% all_present_genes, ]
multiome_sce = multiome_sce[rownames(multiome_sce) %in% all_present_genes, ]


dim(zeb_sce)
dim(all_mouse_sce)
dim(hashikawa_sce_sub )
dim(yalcinbas_sce)
dim(multiome_sce)
#Sanity check
table(rownames(zeb_sce) %in% rownames(all_mouse_sce))
table(rownames(zeb_sce) %in% rownames(hashikawa_sce_sub))
table(rownames(zeb_sce) %in% rownames(yalcinbas_sce))
table(rownames(zeb_sce) %in% rownames(multiome_sce))

#Add the mouse metadata
current_mouse_metadata = readRDS(paste0(wallace_03_data_path, '/wallace_mouse_metaclust_celltype_annot_metadata.rds'))
colData(all_mouse_sce) = S4Vectors::DataFrame(current_mouse_metadata)

full_hashikawa_metadata = readRDS(paste0(hashikawa_data_path, '/hashikawa_mouse_neuron_metaclust_celltype_annot_metadata.rds'))
colData(hashikawa_sce_sub) = full_hashikawa_metadata 

#Rename and add metadata so the columns match across the datasets
#Rename the meta_clust_celltype_annot columns in the zebrafish and mouse data to final_Annotations to match Yalcinbas
zeb_sce$final_Annotations = zeb_sce$meta_clust_celltype_annot
all_mouse_sce$final_Annotations = all_mouse_sce$meta_clust_celltype_annot
hashikawa_sce_sub$final_Annotations = hashikawa_sce_sub$meta_clust_celltype_annot
multiome_sce$final_Annotations = multiome_sce$mid_cluster


zeb_sce$species = 'Zebrafish'
all_mouse_sce$species = 'Wallace|Mouse'
hashikawa_sce_sub$species = 'Hashikawa|Mouse'
yalcinbas_sce$species = 'yalcinbas|Human'
multiome_sce$species = 'multiome|Human'


#Split the mouse and zebrafish SCE objects into lists of SCE objects for each donor, for MetaNeighbor
studies <- unique(all_mouse_sce$study_id)
mouse_sce_list <- lapply(studies, function(study) {
  all_mouse_sce[, all_mouse_sce$study_id == study]
})
names(mouse_sce_list) <- studies


studies <- unique(zeb_sce$study_id)
zebF_sce_list <- lapply(studies, function(study) {
  zeb_sce[, zeb_sce$study_id == study]
})
names(zebF_sce_list) <- studies

#Split Hashikawa by stim or control
studies <- unique(hashikawa_sce_sub$stim)
hashikawa_sce_sub_list <- lapply(studies, function(study) {
  hashikawa_sce_sub[, hashikawa_sce_sub$stim == study]
})
names(hashikawa_sce_sub_list) <- studies


#Get single SCE object
all_donor_sce = mergeSCE(c(mouse_sce_list, zebF_sce_list,hashikawa_sce_sub_list, 
  list(`yalcinbas|Human` = yalcinbas_sce), list(`multiome|Human` = multiome_sce)))
table(all_donor_sce$study_id)
table(all_donor_sce$species)
table(all_donor_sce$final_Annotations)

#Ignore the outlier cells
all_donor_sce = all_donor_sce[, all_donor_sce$final_Annotations!= 'outliers']


#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(
  dat = all_donor_sce,
  min_recurrence = 2,
  exp_labels = all_donor_sce$species
)
length(global_hvgs)
keep_global_hvgs = global_hvgs[1:2000]


MN_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$species,
  cell_type = all_donor_sce$final_Annotations,
  fast_version = TRUE
)

#Plot allby-all AUROC heatmap
plotHeatmap(
  MN_aurocs,
  cex = .5
)
title("MetaNeighbor Human, Mouse, ZebF Habenula: 2000 HVGs")

pdf(paste0(plot_path, '/AllvsAll_MN_human_mouse_zeb_all_cells.pdf'), width = 10, height = 8)
plotHeatmap(
  MN_aurocs,
  cex = .5
)
title("MetaNeighbor Human, Mouse, ZebF Habenula: 2000 HVGs")
dev.off()


#And the best versus next approach

MN_best_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$species,
  cell_type = all_donor_sce$final_Annotations,
  fast_version = TRUE,
  one_vs_best = TRUE,
  symmetric_output = FALSE
)

#Plot best_vs_next AUROC heatmap
plotHeatmap(
  MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Human, Mouse, ZebF Habenula: 2000 HVGs")

pdf(paste0(plot_path, '/BvsNext_MN_human_mouse_zeb_all_cells.pdf'), width = 10, height = 8)
plotHeatmap(
  MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Human, Mouse, ZebF Habenula: 2000 HVGs")
dev.off()


cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)

cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .5)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)

pdf(paste0(plot_path, '/cluster_graph_high_MN_human_mouse_zeb_all_cells.pdf'), width = 10, height = 8)
cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)
dev.off()

pdf(paste0(plot_path, '/cluster_graph_low_MN_human_mouse_zeb_all_cells.pdf'), width = 10, height = 8)
cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .5)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)
dev.off()


#########################
#Just the neurons
##########################


#Just focus on the neuronal populations
all_mouse_sce = all_mouse_sce[, !all_mouse_sce$final_Annotations%in% c('Astrocytes','Polydendrocytes','Differentiating Oligodendrocytes','Oligodendrocytes',
'Pericytes','Fibroblasts','Endothelial','Macrophages', 'Microglia')]
table(all_mouse_sce$final_Annotations)

zeb_sce = zeb_sce[, !zeb_sce$final_Annotations %in% c('non_neuronal')]
table(zeb_sce$final_Annotations)

yalcinbas_sce = yalcinbas_sce[, !yalcinbas_sce$final_Annotations %in% c('Astrocyte', 'Endo','Microglia','Oligo', 'OPC', 'Excit.Thal','Inhib.Thal')]
table(yalcinbas_sce$final_Annotations)

multiome_sce = multiome_sce[, !multiome_sce$final_Annotations %in% c('Astrocyte', 'Endo','Microglia','Oligo', 'OPC', 'Excit.Thal','Inhib.Thal', 'Thal')]
table(multiome_sce$final_Annotations)


#Split the mouse and zebrafish SCE objects into lists of SCE objects for each donor, for MetaNeighbor
studies <- unique(all_mouse_sce$study_id)
mouse_sce_list <- lapply(studies, function(study) {
  all_mouse_sce[, all_mouse_sce$study_id == study]
})
names(mouse_sce_list) <- studies


studies <- unique(zeb_sce$study_id)
zebF_sce_list <- lapply(studies, function(study) {
  zeb_sce[, zeb_sce$study_id == study]
})
names(zebF_sce_list) <- studies

#Split Hashikawa by stim or control
studies <- unique(hashikawa_sce_sub$stim)
hashikawa_sce_sub_list <- lapply(studies, function(study) {
  hashikawa_sce_sub[, hashikawa_sce_sub$stim == study]
})
names(hashikawa_sce_sub_list) <- studies


#Get single SCE object
all_donor_sce = mergeSCE(c(mouse_sce_list, zebF_sce_list,hashikawa_sce_sub_list, 
  list(`yalcinbas|Human` = yalcinbas_sce), list(`multiome|Human` = multiome_sce)))
table(all_donor_sce$study_id)
table(all_donor_sce$species)
table(all_donor_sce$final_Annotations)

#Ignore the outlier cells
all_donor_sce = all_donor_sce[, all_donor_sce$final_Annotations!= 'outliers']


#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(
  dat = all_donor_sce,
  min_recurrence = 2,
  exp_labels = all_donor_sce$species
)
length(global_hvgs)
keep_global_hvgs = global_hvgs[1:2000]


MN_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$species,
  cell_type = all_donor_sce$final_Annotations,
  fast_version = TRUE
)

#Plot allby-all AUROC heatmap
plotHeatmap(
  MN_aurocs,
  cex = .5
)
title("MetaNeighbor Human, Mouse, ZebF Habenula: 2000 HVGs")

pdf(paste0(plot_path, '/AllvsAll_MN_human_mouse_zeb_neurons.pdf'), width = 10, height = 8)
plotHeatmap(
  MN_aurocs,
  cex = .5
)
title("MetaNeighbor Human, Mouse, ZebF Habenula: 2000 HVGs")
dev.off()


#And the best versus next approach

MN_best_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$species,
  cell_type = all_donor_sce$final_Annotations,
  fast_version = TRUE,
  one_vs_best = TRUE,
  symmetric_output = FALSE
)

#Plot best_vs_next AUROC heatmap
plotHeatmap(
  MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Human, Mouse, ZebF Habenula: 2000 HVGs")

pdf(paste0(plot_path, '/BvsNext_MN_human_mouse_zeb_neurons.pdf'), width = 10, height = 8)
plotHeatmap(
  MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Human, Mouse, ZebF Habenula: 2000 HVGs")
dev.off()


cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)

cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .5)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)

pdf(paste0(plot_path, '/cluster_graph_high_MN_human_mouse_zeb_neurons.pdf'), width = 10, height = 8)
cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)
dev.off()

pdf(paste0(plot_path, '/cluster_graph_low_MN_human_mouse_zeb_neurons.pdf'), width = 10, height = 8)
cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .5)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)
dev.off()


#And save the cross-species auroc matrices
saveRDS(MN_aurocs, file = paste0(new_data_path, '/cross_species_MN_AllVsAll_neurons.rds'))

saveRDS(MN_best_aurocs, file = paste0(new_data_path, '/cross_species_MN_BestVsNext_neurons.rds'))

#And save the merged SCE object
saveRDS(all_donor_sce, file = paste0(new_data_path, '/cross_species_merged_sce_neurons.rds'))


#Potential annotations of the multiome clusters

#MHb.1 - Substance P
#MHb1.2 - Cholinergic
#MHb.2 - Cholinergic
#MHb.3 - matches the MHb.3 in Yalcinbas, but not clearly cholinergic or Substance P, and low VGLUT1 compared to the other MHb clusters
#Also doesn't match to the SubP-Cholinergic population in mouse.

#Potentially 3 lateral habenula groups
#The collection that maps to the zebrafish ventral clusters, supposedly the deeply conserved lateral domain
#LHb.1, LHb.1.3.4, LHb.1.3

#The collection that maps to the zebrafish GABA/Glut cluster
#LHb.4

#And the collection that maps to the second mouse LHb cluster
#LHb.2.7


#multiome marker profile

p_bubble = get_bubble_plot_sce(multiome_sce, 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1','TAC3', 'GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1'),
 sample_name = "Multiome Habenula", group_col = "final_Annotations")

p_bubble[[1]]
p_bubble[[2]]

pdf(paste0(plot_path, '/multiome_bubble_Hab_marker_mean.pdf'), width = 10, height = 8)
p_bubble[[1]]
dev.off()

pdf(paste0(plot_path, '/multiome_bubble_Hab_marker_zscore_mean.pdf'), width = 10, height = 8)
p_bubble[[2]]
dev.off()


#And now try just between the human datasets

#Get single SCE object
all_donor_sce = mergeSCE(c(list(`yalcinbas|Human` = yalcinbas_sce), list(`multiome|Human` = multiome_sce)))
table(all_donor_sce$study_id)
table(all_donor_sce$species)
table(all_donor_sce$final_Annotations)

#Ignore the outlier cells
all_donor_sce = all_donor_sce[, all_donor_sce$final_Annotations!= 'outliers']


#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(
  dat = all_donor_sce,
  min_recurrence = 2,
  exp_labels = all_donor_sce$species
)
length(global_hvgs)
keep_global_hvgs = global_hvgs


MN_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$species,
  cell_type = all_donor_sce$final_Annotations,
  fast_version = TRUE
)

#Plot allby-all AUROC heatmap
plotHeatmap(
  MN_aurocs,
  cex = .5
)
title("MetaNeighbor Human Habenula: 894 HVGs")

pdf(paste0(plot_path, '/AllvsAll_MN_human_neurons.pdf'), width = 10, height = 8)
plotHeatmap(
  MN_aurocs,
  cex = .5
)
title("MetaNeighbor Human Habenula: 894 HVGs")
dev.off()


#And the best versus next approach

MN_best_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$species,
  cell_type = all_donor_sce$final_Annotations,
  fast_version = TRUE,
  one_vs_best = TRUE,
  symmetric_output = FALSE
)

#Plot best_vs_next AUROC heatmap
plotHeatmap(
  MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Human Habenula: 894 HVGs")

pdf(paste0(plot_path, '/BvsNext_MN_human_neurons.pdf'), width = 10, height = 8)
plotHeatmap(
  MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Human Habenula: 894 HVGs")
dev.off()


cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)

cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .5)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)

pdf(paste0(plot_path, '/cluster_graph_high_MN_human_neurons.pdf'), width = 10, height = 8)
cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)
dev.off()

pdf(paste0(plot_path, '/cluster_graph_low_MN_human_neurons.pdf'), width = 10, height = 8)
cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .5)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)
dev.off()





#And split by donor
studies <- unique(multiome_sce$orig.ident)
multiome_sce_list <- lapply(studies, function(study) {
  multiome_sce[, multiome_sce$orig.ident == study]
})
names(multiome_sce_list) <- studies

studies <- unique(yalcinbas_sce$Sample)
yalcinbas_sce_list <- lapply(studies, function(study) {
  yalcinbas_sce[, yalcinbas_sce$Sample == study]
})
names(yalcinbas_sce_list) <- studies


#Get single SCE object
all_donor_sce = mergeSCE(c(multiome_sce_list, yalcinbas_sce_list))
table(all_donor_sce$study_id)
table(all_donor_sce$species)
table(all_donor_sce$final_Annotations)

#Ignore the outlier cells
all_donor_sce = all_donor_sce[, all_donor_sce$final_Annotations!= 'outliers']


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
  cell_type = all_donor_sce$final_Annotations,
  fast_version = TRUE
)

#Plot allby-all AUROC heatmap
plotHeatmap(
  MN_aurocs,
  cex = .5
)
title("MetaNeighbor Human cross-donor Habenula: 2000 HVGs")

pdf(paste0(plot_path, '/AllvsAll_MN_human_neurons_cross_donor.pdf'), width = 15, height = 12)
plotHeatmap(
  MN_aurocs,
  cex = .3
)
title("MetaNeighbor Human cross-donor Habenula: 2000 HVGs")
dev.off()


#And the best versus next approach

MN_best_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$final_Annotations,
  fast_version = TRUE,
  one_vs_best = TRUE,
  symmetric_output = FALSE
)

#Plot best_vs_next AUROC heatmap
plotHeatmap(
  MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Human cross-donor Habenula: 2000 HVGs")

pdf(paste0(plot_path, '/BvsNext_MN_human_neurons_cross_donor.pdf'), width = 15, height = 12)
plotHeatmap(
  MN_best_aurocs,
  cex = .3
)
title("MetaNeighbor BvsNext Human cross-donor Habenula: 2000 HVGs")
dev.off()


cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$study_id, all_donor_sce$final_Annotations, size_factor = 3)

cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .5)
plotClusterGraph(cluster_graph, all_donor_sce$study_id, all_donor_sce$final_Annotations, size_factor = 3)

pdf(paste0(plot_path, '/cluster_graph_high_MN_human_neurons_cross_donor.pdf'), width = 10, height = 8)
cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$study_id, all_donor_sce$final_Annotations, size_factor = 3)
dev.off()

pdf(paste0(plot_path, '/cluster_graph_low_MN_human_neurons_cross_donor.pdf'), width = 10, height = 8)
cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .5)
plotClusterGraph(cluster_graph, all_donor_sce$study_id, all_donor_sce$final_Annotations, size_factor = 3)
dev.off()



#The MHb3 cluster is distinct and replicable across the human datasets, a bit weak though, since the Yalcinbas cluster is only like 15 cells
#But distinct across donors in the multiome, a bit difficult there since the clusters were defined included all donors, but still seems strong


#So all in all, what I would label clusters
#MHb Substance P
#MHb Cholinergic
#MHb undefined. This is the MHb3 in humans. And decidely does not map to the SubP_Cholinergic cluster in mouse

#LHb1 - the Lateral Habenula clusters that map to the zebrafish ventral clusters
#LHb2 - the lateral Habenula clusters that map to the zebrafish GABA/Glut cluster

#LHb1 and 2 seem related, but the GABA/Glut clusters form the strongest cross-species cluster

#LHb3 - a potential mammal-specific lateral cluster, Maps to the Wallace LHb2, Yalcinbas 2 and 7, and multiome, 2.7 and 7


#Get a barplot of cell numbers per cluster, than a stacked barplot of donor contribution to each cluster

cell_num_df = as.data.frame(colData(multiome_sce)) %>% select(orig.ident, final_Annotations) %>% group_by(final_Annotations, orig.ident) %>% 
  summarise(n = n())

annotation_order <- cell_num_df %>%
  group_by(final_Annotations) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(final_Annotations)

# Stacked barplot with proportions
p1 = cell_num_df %>%
  group_by(final_Annotations) %>%
  mutate(prop = n / sum(n),
         final_Annotations = factor(final_Annotations, levels = annotation_order)) %>%
  ggplot(aes(x = final_Annotations, y = prop, fill = orig.ident)) +
  geom_col() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Habenula clusters", y = "Proportion", fill = "Sample") +
  scale_fill_manual(values = MetBrewer::met.brewer("Redon", n = 10), name = 'Donor')

# Total cell count barplot
p2 = cell_num_df %>%
  group_by(final_Annotations) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  mutate(final_Annotations = factor(final_Annotations, levels = annotation_order)) %>%
  ggplot(aes(x = final_Annotations, y = total, fill = final_Annotations)) +
  geom_col() +
  ggbreak::scale_y_break(c(5000, 10000)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Habenula Clusters", y = "Total Cell Count") +
  scale_fill_manual(values = MetBrewer::met.brewer("Navajo", n = 10), name = 'Habenula clusters')


ggsave(plot = p1, filename = 'cluster_per_donor_stacked_bar.pdf', path = plot_path, 
width = 4, height = 4,device = 'pdf', useDingbats = F)

ggsave(plot = p2, filename = 'cluster_total_cell_number_bar.pdf', path = plot_path, 
width = 4, height = 4,device = 'pdf', useDingbats = F, onefile = FALSE)




#Check out some DE and markers
#Check out the markers when including all the data
all_multiome_markers = compute_markers(assay(multiome_sce, 'cpm'), multiome_sce$mid_cluster)
all_multiome_markers %>% filter(gene %in% c('ESR1', 'PVALB', 'KIT')) %>% group_by(gene) %>% arrange(average_expression, .by_group = T) %>% View()

table(multiome_sce$mid_cluster, multiome_sce$orig.ident)
table(multiome_sce$mid_cluster)


p_bubble = get_bubble_plot_sce(multiome_sce, 
  top_markers = c('ESR1', 'PVALB','KIT', 'GAD1', 'GAD2', 'SLC32A1', 'SLC17A7', 'OPRM1'),
 sample_name = "Multiome Habenula", group_col = "final_Annotations")

p_bubble[[1]]
p_bubble[[2]]



ggplot(all_multiome_markers %>% filter(cell_type == 'MHb.1'), aes(x = log2(fold_change), y = auroc)) + 
  geom_point(color = 'grey', alpha = .5) +
  theme_bw() + ggtitle('Multiome MHb.1 DE') + ylim(0,1) + xlim(-7.5, 5.5) +
  geom_point(data = all_multiome_markers %>% filter(cell_type == 'MHb.1' & gene %in% c('TAC1', 'TAC3', 'CHAT', 'SLC5A7', 'SLC18A3')),
             color = "red", size = 3) +
  geom_label_repel(data = all_multiome_markers %>% filter(cell_type == 'MHb.1' & gene %in% c('TAC1', 'TAC3', 'CHAT', 'SLC5A7', 'SLC18A3')),
                   aes(label = gene),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.color = "black",
                   segment.size = 0.5)


ggplot(all_multiome_markers %>% filter(cell_type == 'MHb.2'), aes(x = log2(fold_change), y = auroc)) + 
  geom_point(color = 'grey', alpha = .5) +
  theme_bw() + ggtitle('Multiome MHb.2 DE') + ylim(0,1) + xlim(-7.5, 5.5) +
  geom_point(data = all_multiome_markers %>% filter(cell_type == 'MHb.2' & gene %in% c('TAC1', 'TAC3', 'CHAT', 'SLC5A7', 'SLC18A3')),
             color = "red", size = 3) +
  geom_label_repel(data = all_multiome_markers %>% filter(cell_type == 'MHb.2' & gene %in% c('TAC1', 'TAC3', 'CHAT', 'SLC5A7', 'SLC18A3')),
                   aes(label = gene),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.color = "black",
                   segment.size = 0.5)

ggplot(all_multiome_markers %>% filter(cell_type == 'MHb.3'), aes(x = log2(fold_change), y = auroc)) + 
  geom_point(color = 'grey', alpha = .5) +
  theme_bw() + ggtitle('Multiome MHb.3 DE')  + ylim(0,1) + xlim(-7.5, 5.5) +
  geom_point(data = all_multiome_markers %>% filter(cell_type == 'MHb.3' & gene %in% c('TAC1', 'TAC3', 'CHAT', 'SLC5A7', 'SLC18A3')),
             color = "red", size = 3) +
  geom_label_repel(data = all_multiome_markers %>% filter(cell_type == 'MHb.3' & gene %in% c('TAC1', 'TAC3', 'CHAT', 'SLC5A7', 'SLC18A3')),
                   aes(label = gene),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.color = "black",
                   segment.size = 0.5)


ggplot(all_multiome_markers %>% filter(cell_type == 'LHb.4'), aes(x = log10(average_expression), y = auroc)) + 
  geom_point(color = 'grey', alpha = .5) +
  theme_bw() + ggtitle('Multiome LHb.4 DE') +
  geom_point(data = all_multiome_markers %>% filter(cell_type == 'LHb.4' & gene %in% c('GAD1', 'GAD2', 'SLC32A1', 'GPR151','ESR1', 'PVALB', 'KIT', 'SLC17A6', 'SLC17A7')),
             color = "red", size = 3) +
  geom_label_repel(data = all_multiome_markers %>% filter(cell_type == 'LHb.4' & gene %in% c('GAD1', 'GAD2', 'SLC32A1', 'GPR151','ESR1', 'PVALB', 'KIT', 'SLC17A6', 'SLC17A7')),
                   aes(label = gene),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.color = "black",
                   segment.size = 0.5)


ggplot(all_multiome_markers %>% filter(cell_type == 'LHb.2.7'), aes(x = log2(fold_change), y = auroc)) + 
  geom_point(color = 'grey', alpha = .5) +
  theme_bw() + ggtitle('Multiome LHb.2.7 DE') +
  geom_point(data = all_multiome_markers %>% filter(cell_type == 'LHb.2.7' & gene %in% c('GAD1', 'GAD2', 'SLC32A1', 'GPR151','ESR1', 'PVALB', 'KIT', 'SLC17A6', 'SLC17A7', 'OPRM1')),
             color = "red", size = 3) +
  geom_label_repel(data = all_multiome_markers %>% filter(cell_type == 'LHb.2.7' & gene %in% c('GAD1', 'GAD2', 'SLC32A1', 'GPR151','ESR1', 'PVALB', 'KIT', 'SLC17A6', 'SLC17A7', 'OPRM1')),
                   aes(label = gene),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.color = "black",
                   segment.size = 0.5)

ggplot(all_multiome_markers %>% filter(cell_type == 'LHb.1'), aes(x = log2(fold_change), y = auroc)) + 
  geom_point(color = 'grey', alpha = .5) +
  theme_bw() + ggtitle('Multiome LHb.1 DE') +
  geom_point(data = all_multiome_markers %>% filter(cell_type == 'LHb.1' & gene %in% c('GAD1', 'GAD2', 'SLC32A1', 'GPR151','GAP43', 'SLC17A6', 'SLC17A7')),
             color = "red", size = 3) +
  geom_label_repel(data = all_multiome_markers %>% filter(cell_type == 'LHb.1' & gene %in% c('GAD1', 'GAD2', 'SLC32A1', 'GPR151','GAP43', 'SLC17A6', 'SLC17A7')),
                   aes(label = gene),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   segment.color = "black",
                   segment.size = 0.5)








