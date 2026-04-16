#Rerunning the same cross-species assessment, but with the refined mid resolution multiome clusters

library(SingleCellExperiment)
library(MetaNeighbor)
library(qs2)
library(dplyr)
library(ggplot2)
library(here)


here::here()

#Path for new data generated
new_data_path = here('processed-data', '05_03_annotation_adjustments','07_redo_cross_species')
#Path to plot directory
plot_path = here('plots', '05_03_annotation_adjustments','07_redo_cross_species')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)

#This is the same multiome SCE as before, but with the refined annotations added in
#in refined_mid_cluster metadata
multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

multiome_sce
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))


#Path to the Yalcinbas pilot data
#Going with the official_final_sce.RDATA
yalcinbas_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/sce_objects'

#Nonhuman data paths
zeb_data_path = here('processed-data', '05_02_external_Hb_comparisons','07_cross_species_hab')
mouse_data_path = here('processed-data', '05_02_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')
wallace_03_data_path = here('processed-data', '05_02_external_Hb_comparisons', '03_metaMarkers_wallace_2019')

hashikawa_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data'
hashikawa_data_path = here('processed-data','05_02_external_Hb_comparisons','04_wallace_hashikawa_mouse')

#Path to orthologs
path_to_orthologs = here('processed-data', '05_02_external_Hb_comparisons', 'human_mouse_zebrafish_orthologs.txt.gz')



#Yalcinbas data
#Loads as an object labeled 'sce'
#16437 cells
load(paste0(yalcinbas_path, '/official_final_sce.RDATA'))
yalcinbas_sce = sce
rm(sce)

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

#Clean memory
gc()

# Add the mouse metadata
current_mouse_metadata = readRDS(paste0(wallace_03_data_path, '/wallace_mouse_metaclust_celltype_annot_metadata.rds'))
colData(all_mouse_sce) = S4Vectors::DataFrame(current_mouse_metadata)

full_hashikawa_metadata = readRDS(paste0(hashikawa_data_path, '/hashikawa_mouse_neuron_metaclust_celltype_annot_metadata.rds'))
colData(hashikawa_sce_sub) = full_hashikawa_metadata 

#Rename and add metadata so the columns match across the datasets
#Rename the meta_clust_celltype_annot columns in the zebrafish and mouse data to final_Annotations to match Yalcinbas
zeb_sce$final_Annotations = zeb_sce$meta_clust_celltype_annot
all_mouse_sce$final_Annotations = all_mouse_sce$meta_clust_celltype_annot
hashikawa_sce_sub$final_Annotations = hashikawa_sce_sub$meta_clust_celltype_annot
multiome_sce$final_Annotations = multiome_sce$refined_mid_cluster


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

multiome_sce = multiome_sce[, !multiome_sce$final_Annotations %in% c('Astrocyte', 'Endo','Microglia','Oligo', 'OPC', 'Excit.Thal','Inhib.Thal', 'Thal', 'Ependymal')]
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



#Make the default cluster graph a bit nicer

cluster_graph = makeClusterGraph(MN_best_aurocs , low_threshold = .5)

plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)

#Make the default graph a bit nicer
#Species by color
#Dataset gets its own shape

metadata_df = as.data.frame(colData(all_donor_sce))
metadata_df$species <- gsub("\\|", ".", metadata_df$species)
metadata_df %>% View()


temp_cluster_graph = cluster_graph

species_color_palette = c('Wallace.Mouse' = 'brown','Hashikawa.Mouse' = 'brown',
 'yalcinbas.Human' = 'purple', 'multiome.Human' = 'purple', 'Zebrafish' = 'darkblue')

species_shape_palette = c('Wallace.Mouse' = 'square',
                          'Hashikawa.Mouse' = 'circle',
                          'yalcinbas.Human' = 'square',
                          'multiome.Human' = 'circle',
                          'Zebrafish' = 'circle')

#Subclass colors for the mouse clusters
species_colors = species_color_palette[metadata_df$species]
names(species_colors) = paste(metadata_df$species, metadata_df$final_Annotations, sep = '|')

vertex_study <- factor(igraph::V(temp_cluster_graph)$name)

#Set up the vertex colors
levels(vertex_study) <- species_colors[levels(vertex_study)]

igraph::V(temp_cluster_graph)$color <- as.character(vertex_study)
igraph::V(temp_cluster_graph)$label.color <- "black"


# Create shape vector similar to species_colors
species_shapes = species_shape_palette[metadata_df$species]
names(species_shapes) = paste(metadata_df$species, metadata_df$final_Annotations, sep = '|')

# Set up the vertex shapes
vertex_study_shape <- factor(igraph::V(temp_cluster_graph)$name)
levels(vertex_study_shape) <- species_shapes[levels(vertex_study_shape)]

igraph::V(temp_cluster_graph)$shape <- as.character(vertex_study_shape)


#Just the cluster name as the label, excluding the study
full_vertex_labels = igraph::V(temp_cluster_graph)$name
vertex_labels = sapply(strsplit(igraph::V(temp_cluster_graph)$name, split = '|', fixed = T), '[[', 2)
igraph::V(temp_cluster_graph)$label <- vertex_labels


igraph::E(temp_cluster_graph)$width <- igraph::E(temp_cluster_graph)$weight
igraph::E(temp_cluster_graph)$color <- c("orange","darkgray")[as.numeric(igraph::E(temp_cluster_graph)$weight >= 0.5) + 1]

size_factor=1
plot(temp_cluster_graph,
     vertex.size = 8,
     vertex.label.cex=1.5,
     vertex.label.family="Helvetica",
     vertex.label.font=1,
     vertex.frame.color = 'black',
     edge.width = igraph::E(temp_cluster_graph)$width * 5,
     edge.arrow.size=.1*5,
     edge.arrow.width=0.5*5)
graphics::legend(
    "topleft", 
    legend = names(species_color_palette), 
    pt.bg = species_color_palette,
    col = species_color_palette,
    pt.cex = 2, cex = 0.5*2, pch=c(15, 16, 15, 16, 16)
)

pdf(paste0(plot_path, '/clean_cluster_graph_low_MN_human_mouse_zeb_neurons.pdf'), width = 10, height = 8)
plot(temp_cluster_graph,
     vertex.size = 8,
     vertex.label.cex=1.5,
     vertex.label.family="Helvetica",
     vertex.label.font=1,
     vertex.frame.color = 'black',
     edge.width = igraph::E(temp_cluster_graph)$width * 5,
     edge.arrow.size=.1*5,
     edge.arrow.width=0.5*5)
graphics::legend(
    "topleft", 
    legend = names(species_color_palette), 
    pt.bg = species_color_palette,
    col = species_color_palette,
    pt.cex = 2, cex = 0.5*2, pch=c(15, 16, 15, 16, 16)
)
dev.off()





