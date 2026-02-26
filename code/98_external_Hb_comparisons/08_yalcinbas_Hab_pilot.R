#Loading up and seeing what is currently present in the published Yalcinbas Habenula single cell data
#The single cell data is all neurotypical donors


library(SingleCellExperiment)
library(Seurat)
library(MetaNeighbor)
library(dplyr)
library(ggplot2)
library(here)

here::here()


#Path to the Yalcinbas pilot data
#Going with the official_final_sce.RDATA
yalcinbas_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/sce_objects'
list.files(yalcinbas_path)

#Nonhuman data paths
zeb_data_path = here('processed-data', '98_external_Hb_comparisons','07_cross_species_hab')
mouse_data_path = here('processed-data', '98_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')
wallace_03_data_path = here('processed-data', '98_external_Hb_comparisons', '03_metaMarkers_wallace_2019')

hashikawa_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data'
hashikawa_data_path = here('processed-data','98_external_Hb_comparisons','04_wallace_hashikawa_mouse')

#Path to orthologs
path_to_orthologs = here('processed-data', '98_external_Hb_comparisons', 'human_mouse_zebrafish_orthologs.txt.gz')


#Path for new data generated
new_data_path = here('processed-data', '98_external_Hb_comparisons','08_yalcinbas_Hab_pilot')
#Path to plot directory
plot_path = here('plots', '98_external_Hb_comparisons', '08_yalcinbas_Hab_pilot')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)

#Source the bubble plot functions
source(here('code','98_external_Hb_comparisons', 'bubble_plot_functions.R'))




#Loads as an object labeled 'sce'
#16437 cells
load(paste0(yalcinbas_path, '/official_final_sce.RDATA'))
yalcinbas_sce = sce
rm(sce)

colnames(colData(yalcinbas_sce))
table(yalcinbas_sce$final_Annotations)
table(yalcinbas_sce$broad_Annotations)

table(yalcinbas_sce$Sample, yalcinbas_sce$final_Annotations )

yalcinbas_sce

assay(yalcinbas_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(yalcinbas_sce, 'counts'))


#Just try the whole dataset as a sample, can break it up by donor later, but there are only 4 donors with decent annotated habenula cells

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
all_present_genes = Reduce(intersect, list(rownames(all_mouse_sce), rownames(hashikawa_sce_sub), rownames(zeb_sce), rownames(yalcinbas_sce)))
all_present_genes = all_present_genes[!is.na(all_present_genes)]
length(all_present_genes)


#Filter the SCE objects for the present genes
zeb_sce = zeb_sce[rownames(zeb_sce) %in% all_present_genes, ]
all_mouse_sce = all_mouse_sce[rownames(all_mouse_sce) %in% all_present_genes, ]
hashikawa_sce_sub = hashikawa_sce_sub[rownames(hashikawa_sce_sub) %in% all_present_genes, ]
yalcinbas_sce = yalcinbas_sce[rownames(yalcinbas_sce) %in% all_present_genes, ]

dim(zeb_sce)
dim(all_mouse_sce)
dim(hashikawa_sce_sub )
dim(yalcinbas_sce)
#Sanity check
table(rownames(zeb_sce) %in% rownames(all_mouse_sce))
table(rownames(zeb_sce) %in% rownames(hashikawa_sce_sub))
table(rownames(zeb_sce) %in% rownames(yalcinbas_sce))


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

zeb_sce$species = 'Zebrafish'
all_mouse_sce$species = 'Wallace|Mouse'
hashikawa_sce_sub$species = 'Hashikawa|Mouse'
yalcinbas_sce$species = 'Human'


#Just focus on the neuronal populations
all_mouse_sce = all_mouse_sce[, !all_mouse_sce$final_Annotations%in% c('Astrocytes','Polydendrocytes','Differentiating Oligodendrocytes','Oligodendrocytes',
'Pericytes','Fibroblasts','Endothelial','Macrophages', 'Microglia')]
table(all_mouse_sce$final_Annotations)

zeb_sce = zeb_sce[, !zeb_sce$final_Annotations %in% c('non_neuronal')]
table(zeb_sce$final_Annotations)

yalcinbas_sce = yalcinbas_sce[, !yalcinbas_sce$final_Annotations %in% c('Astrocyte', 'Endo','Microglia','Oligo', 'OPC', 'Excit.Thal','Inhib.Thal')]
table(yalcinbas_sce$final_Annotations)


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
all_donor_sce = mergeSCE(c(mouse_sce_list, zebF_sce_list,hashikawa_sce_sub_list, list(Human = yalcinbas_sce)))
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




#From all the cross-species, inital annotation of Yalcinbas

#LHb 1, 3, 4, 5 all map to the VGLUT2 + inhibitory populations across species
#LHb 2 and 7 map to the other lateral Hab group
#LHb 6 is kind of off on its own

#MHb2 is cholinergic
#MHb1 is substance P

#MHb3 is off on its own, there are only 18 cells here too



#And what about just mouse and human

#Reload to focus on the human and mouse 1to1 orthologs
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

#Yalcinbas data
load(paste0(yalcinbas_path, '/official_final_sce.RDATA'))
yalcinbas_sce = sce
rm(sce)

assay(yalcinbas_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(yalcinbas_sce, 'counts'))


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

table(rownames(all_mouse_sce) %in% rownames(hashikawa_sce_sub) )

#Get the set of genes that are present in all datasets, and are 1-to-1 orthologs
#12412 across all three, which is pretty good
all_present_genes = Reduce(intersect, list(rownames(all_mouse_sce), rownames(hashikawa_sce_sub), rownames(yalcinbas_sce)))
all_present_genes = all_present_genes[!is.na(all_present_genes)]
length(all_present_genes)


#Filter the SCE objects for the present genes
all_mouse_sce = all_mouse_sce[rownames(all_mouse_sce) %in% all_present_genes, ]
hashikawa_sce_sub = hashikawa_sce_sub[rownames(hashikawa_sce_sub) %in% all_present_genes, ]
yalcinbas_sce = yalcinbas_sce[rownames(yalcinbas_sce) %in% all_present_genes, ]


dim(all_mouse_sce)
dim(hashikawa_sce_sub )
dim(yalcinbas_sce)
#Sanity check
table(rownames(all_mouse_sce) %in% rownames(hashikawa_sce_sub))
table(rownames(all_mouse_sce) %in% rownames(yalcinbas_sce))


#Add the mouse metadata
current_mouse_metadata = readRDS(paste0(wallace_03_data_path, '/wallace_mouse_metaclust_celltype_annot_metadata.rds'))
colData(all_mouse_sce) = S4Vectors::DataFrame(current_mouse_metadata)

full_hashikawa_metadata = readRDS(paste0(hashikawa_data_path, '/hashikawa_mouse_neuron_metaclust_celltype_annot_metadata.rds'))
colData(hashikawa_sce_sub) = full_hashikawa_metadata 

#Rename and add metadata so the columns match across the datasets
all_mouse_sce$final_Annotations = all_mouse_sce$meta_clust_celltype_annot
hashikawa_sce_sub$final_Annotations = hashikawa_sce_sub$meta_clust_celltype_annot

all_mouse_sce$species = 'Wallace|Mouse'
hashikawa_sce_sub$species = 'Hashikawa|Mouse'
yalcinbas_sce$species = 'Human'


#Just focus on the neuronal populations
all_mouse_sce = all_mouse_sce[, !all_mouse_sce$final_Annotations%in% c('Astrocytes','Polydendrocytes','Differentiating Oligodendrocytes','Oligodendrocytes',
'Pericytes','Fibroblasts','Endothelial','Macrophages', 'Microglia')]
table(all_mouse_sce$final_Annotations)

yalcinbas_sce = yalcinbas_sce[, !yalcinbas_sce$final_Annotations %in% c('Astrocyte', 'Endo','Microglia','Oligo', 'OPC', 'Excit.Thal','Inhib.Thal')]
table(yalcinbas_sce$final_Annotations)


#Get single SCE object
human_mouse_sce = mergeSCE(c(mouse_sce_list, hashikawa_sce_sub_list, list(Human = yalcinbas_sce)))
table(human_mouse_sce$study_id)
table(human_mouse_sce$species)
table(human_mouse_sce$final_Annotations)

#Ignore the outlier cells
human_mouse_sce = human_mouse_sce[, human_mouse_sce$final_Annotations!= 'outliers']


#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(
  dat = human_mouse_sce,
  min_recurrence = 2,
  exp_labels = human_mouse_sce$species
)
length(global_hvgs)
#keep_global_hvgs = global_hvgs[1:2000]
keep_global_hvgs = global_hvgs

huMous_MN_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = human_mouse_sce,
  study_id = human_mouse_sce$species,
  cell_type = human_mouse_sce$final_Annotations,
  fast_version = TRUE
)

#Plot allby-all AUROC heatmap
plotHeatmap(
  huMous_MN_aurocs,
  cex = .5
)
title("MetaNeighbor Human, Mouse Habenula: 1890 HVGs")

pdf(paste0(plot_path, '/AllvsAll_MN_human_mouse_neurons.pdf'), width = 10, height = 8)
plotHeatmap(
  huMous_MN_aurocs,
  cex = .5
)
title("MetaNeighbor Human, Mouse Habenula: 1890 HVGs")
dev.off()


#And the best versus next approach

huMous_MN_best_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = human_mouse_sce,
  study_id = human_mouse_sce$species,
  cell_type = human_mouse_sce$final_Annotations,
  fast_version = TRUE,
  one_vs_best = TRUE,
  symmetric_output = FALSE
)

#Plot best_vs_next AUROC heatmap
plotHeatmap(
  huMous_MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Human, Mouse Habenula: 1890 HVGs")


pdf(paste0(plot_path, '/BvsNext_MN_human_mouse_neurons.pdf'), width = 10, height = 8)
plotHeatmap(
  huMous_MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Human, Mouse Habenula: 1890 HVGs")
dev.off()

cluster_graph = makeClusterGraph(huMous_MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, human_mouse_sce$species, human_mouse_sce$final_Annotations, size_factor = 3)


pdf(paste0(plot_path, '/cluster_graph_MN_human_mouse_neurons.pdf'), width = 10, height = 8)
cluster_graph = makeClusterGraph(huMous_MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, human_mouse_sce$species, human_mouse_sce$final_Annotations, size_factor = 3)
dev.off()

#Again, the predictions from zebrafish to human hold up, are a bit more obvious when comparing human and mouse
#LHb 1, 3, 4, 5 all map to the VGLUT2 + inhibitory populations across species
#LHb 2 and 7 map to the other lateral Hab group
#LHb 6 is kind of off on its own

#MHb2 is cholinergic
#MHb1 is substance P

#MHb3 is off on its own, there are only 18 cells here too



#Check out the marker expression in the Yalcinbas SCE object, using the final_Annotations

p_bubble = get_bubble_plot_sce(yalcinbas_sce, 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1','TAC3', 'GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1'),
 sample_name = "Yalcinbas Habenula", group_col = "final_Annotations")

p_bubble[[1]]
p_bubble[[2]]

pdf(paste0(plot_path, '/yalcinbas_bubble_Hab_marker_mean.pdf'), width = 10, height = 8)
p_bubble[[1]]
dev.off()

pdf(paste0(plot_path, '/yalcinbas_bubble_Hab_marker_zscore_mean.pdf'), width = 10, height = 8)
p_bubble[[2]]
dev.off()


#Compute human markers
markers_human = MetaMarkers::compute_markers(assay(yalcinbas_sce, "cpm"), yalcinbas_sce$final_Annotations)

markers_human %>% group_by(cell_type) %>% slice_max(order_by = auroc, n = 25) %>% View()


#Ignoring the iffy clusters

p_bubble = get_bubble_plot_sce(yalcinbas_sce[ , !yalcinbas_sce$final_Annotations %in% c('LHb.6', 'MHb.3')], 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1','TAC3', 'GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1'),
 sample_name = "Yalcinbas Habenula", group_col = "final_Annotations")

p_bubble[[1]]
p_bubble[[2]]



#Alright, so the simple split of LHb into the inhib and nonInhib expressing groups doesnt seem to be as clean in humans.
#To recap, Wallace mouse had 2 Lateral clusters, those were distinguished by + or - GAD1, 2, and VGAT expression.
#Hashikawa also had some LHb GAD1, 2 and VGAT expressing clusters map to the same GAD1, 2, and VGAT in Wallace, but there were also some non LHb inhib expressing cells that mapped
#The Zebrafish GAP43 and GAD2, 1, VGAT cluster also maps here

#In the Yalcinbase dataset, LHb7 and LHb4 have that mixed excitatory and inhibitory signature, but they map to different mouse Lateral groups.
#So maybe the defining feature of the two lateral hab groups is not simply the excit/inhib co-expression


#Still, get a broader grouping of the Yalcinbas by the cross-species replicability and check out their markers.
#Essentially the cholinergic and substance P for medial, and then Lateral 1 and Lateral 2 for lateral.

#Ignore MHb3 and LHb6 for now

meta_annot_vec = c( 'SupstanceP_medial' = 'MHb.1',
                    'Cholinergic_medial' = 'MHb.2',
                    'Lateral_1' = 'LHb.7',
                    'Lateral_1' = 'LHb.2',
                    'Lateral_2' = 'LHb.1',
                    'Lateral_2' = 'LHb.3',
                    'Lateral_2' = 'LHb.4',
                    'Lateral_2' = 'LHb.5',
                    'outliers' = 'MHb.3' ,
                    'outliers' = 'LHb.6'   
)

meta_annot_vec  = setNames(names(meta_annot_vec), meta_annot_vec)


yalcinbas_sce$broad_annot = unname(meta_annot_vec[yalcinbas_sce$final_Annotations])
table(yalcinbas_sce$broad_annot,yalcinbas_sce$final_Annotations)


#Compute human markers
broad_markers_human = MetaMarkers::compute_markers(assay(yalcinbas_sce, "cpm"), yalcinbas_sce$broad_annot)

broad_markers_human %>% group_by(cell_type) %>% slice_max(order_by = fold_change, n = 25) %>% View()








