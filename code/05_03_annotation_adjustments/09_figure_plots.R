#Pulling together plots for figures


library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(qs2)
library(scater)
library(here)

here::here()


#Path to save any generated data
new_data_path = here('processed-data', '05_03_annotation_adjustments', '09_figure_plots')
#Path to plot directory
plot_path = here('plots','05_03_annotation_adjustments', '09_figure_plots')


if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))
colnames(colData(multiome_sce))

#SingleCellExperiment::altExp(multiome_sce, "ATAC") <- NULL

#bubble plots
source(here('code','05_02_external_Hb_comparisons', 'bubble_plot_functions.R'))

#colors
source(here('code','05_03_annotation_adjustments','celltype_colors.R'))




#First are the fine to mid cluster adjustments that should have been done before any cross-species
#This is the C.21.Astrocyte being called ependymal and the C.11.MHb.1.2 being called MHb2 based on cholinergic markers
#Show those targeted bubble plots first


get_bubble_plot_sce(multiome_sce,
  top_markers = c('CHAT','SLC18A3','SLC5A7'),
  sample_name = "Multiome Habenula", group_col = "cluster_ann"
)


cell_filt = multiome_sce$mid_cluster %in% c('MHb.3','MHb.2','MHb.1','MHb.1.2')
celltype_order = c('C.36.MHb.3','C.07.MHb.2','C.11.MHb.1.2','C.16.MHb.1.2','C.10.MHb.1','C.14.MHb.1')

p_bubble = get_bubble_plot_sce(multiome_sce[ , cell_filt], 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3'),
 sample_name = "Multiome Habenula: Medial habenula fine clusters", group_col = "cluster_ann", group_order = celltype_order)

p_bubble[[1]]
p_bubble[[2]] 

ggsave(p_bubble[[2]], filename = paste0(plot_path, '/MHb_cholinergic_bubble_plot_zscore.pdf'), device = 'pdf', 
width = 6, height = 4)


#And now just the non-neuronal populations
cell_filt = multiome_sce$mid_cluster %in% c('Astrocyte','Oligo','OPC','Microglia','Endo')
celltype_order = c('C.21.Astrocyte','C.20.Astrocyte','C.02.Oligo','C.22.Oligo','C.26.OPC',
'C.27.Microglia','C.41.Microglia','C.29.Endo')

p_bubble = get_bubble_plot_sce(multiome_sce[ , cell_filt], 
  top_markers = c('FOXJ1','TPPP3','PIFO','CD24','TMEM212','RARRES2','GFAP','AQP4','SOX9','HOPX','SLC1A2',
  'TJP1','S100B','MOG','SOX10','OLIG1','OLIG2','PDGFRA','P2RY12','TMEM119','CX3CR1','CD68','CSF1R', 'CLDN5', 'CDH5', 'PECAM1'),
 sample_name = "Multiome Habenula: Non-neuronal clusters", group_col = "cluster_ann", group_order = celltype_order)

p_bubble[[1]]
p_bubble[[2]] 

ggsave(p_bubble[[2]], filename = paste0(plot_path, '/non_neurons_bubble_plot_zscore.pdf'), device = 'pdf', 
width = 8, height = 4)


#Currently, those adjustments are just reflected in the refined_mid_cluster annotations
#Add them to the mid_cluster annotations too for the figure
multiome_sce$mid_cluster_adj = multiome_sce$mid_cluster
multiome_sce$mid_cluster_adj[multiome_sce$cluster_ann == 'C.21.Astrocyte'] = 'Ependymal'
multiome_sce$mid_cluster_adj[multiome_sce$cluster_ann == 'C.11.MHb.1.2'] = 'MHb.2'


#General umaps with the annotations before any cross-species comparisons

celltype_order = c('MHb.1','MHb.1.2','MHb.2', 'MHb.3',
'LHb.2.7','LHb.1','LHb.1.3','LHb.1.3.4','LHb.4',
'Inhib.Thal','Excit.Thal',
'Ependymal','Astrocyte','Microglia','Endo','Oligo','OPC')

cons_celltype_order = c('MHb_A','MHb_C','MHb_B', 'MHb_D',
'LHb_A','LHb_B','LHb_C', 'GABA_LHb_C.1', 'GABA_LHb_C.2', 
'Inhib.Thal','Excit.Thal',
'Ependymal','Astrocyte','Microglia','Endo','Oligo','OPC')

multiome_sce$mid_cluster_adj <- factor(multiome_sce$mid_cluster_adj, levels = celltype_order)
multiome_sce$refined_mid_cluster <- factor(multiome_sce$refined_mid_cluster, levels = cons_celltype_order)

plt1 <- plotReducedDim(multiome_sce, 
                dimred = "wnn.umap",
                colour_by = "mid_cluster_adj", 
              point_size = .5) +
  scale_color_manual(values = my_colors_mid) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), reverse = TRUE)) +
  labs(
    title = "WNN (RNA+ATAC) UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Cell type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

plt1

plt2 <- plotReducedDim(multiome_sce, 
                dimred = "umap.integrated",
                colour_by = "mid_cluster_adj", 
              point_size = .5) +
  scale_color_manual(values = my_colors_mid) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), reverse = TRUE)) +
  labs(
    title = "RNA UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Cell type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

plt2

plt3 <- plotReducedDim(multiome_sce, 
                dimred = "umap.lsi.integrated",
                colour_by = "mid_cluster_adj", 
              point_size = .5) +
  scale_color_manual(values = my_colors_mid) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), reverse = TRUE)) +
  labs(
    title = "ATAC UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Cell type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

plt3

#And then the rna UMAP with the final annotations

plt4 <- plotReducedDim(multiome_sce, 
                dimred = "umap.integrated",
                colour_by = "refined_mid_cluster", 
              point_size = .5) +
  scale_color_manual(values = my_colors_mid) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), reverse = TRUE)) +
  labs(
    title = "RNA UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Cell type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )
# Build centroid table for labels
# emb <- as.data.frame(reducedDim(multiome_sce, "umap.integrated"))
# colnames(emb)[1:2] <- c("UMAP1", "UMAP2")
# emb$cluster <- multiome_sce$refined_mid_cluster

# centers <- emb |>
#   summarize(
#     UMAP1 = median(UMAP1),
#     UMAP2 = median(UMAP2),
#     .by = cluster
#   )

# # Add labels
# plt4 = plt4 +
#   ggrepel::geom_text_repel(
#     data = centers,
#     aes(x = UMAP1, y = UMAP2, label = cluster),
#     inherit.aes = FALSE,
#     size = 3
#   )


plt4


ggsave(plt1, filename = 'wnn_umap_starting_mid_cluster.pdf', path = plot_path, device = 'pdf', 
height = 6, width = 6)
ggsave(plt2, filename = 'RNA_umap_starting_mid_cluster.pdf', path = plot_path, device = 'pdf', 
height = 6, width = 6)
ggsave(plt3, filename = 'ATAC_umap_starting_mid_cluster.pdf', path = plot_path, device = 'pdf', 
height = 6, width = 6)
ggsave(plt4, filename = 'RNA_umap_final_mid_cluster.pdf', path = plot_path, device = 'pdf', 
height = 6, width = 6)
#pngs
ggsave(plt1, filename = 'wnn_umap_starting_mid_cluster.png', path = plot_path, device = 'png', 
height = 6, width = 6)
ggsave(plt2, filename = 'RNA_umap_starting_mid_cluster.png', path = plot_path, device = 'png', 
height = 6, width = 6)
ggsave(plt3, filename = 'ATAC_umap_starting_mid_cluster.png', path = plot_path, device = 'png', 
height = 6, width = 6)
ggsave(plt4, filename = 'RNA_umap_final_mid_cluster.png', path = plot_path, device = 'png', 
height = 6, width = 6)

#And pull together the broader marker bubble plot, here we want the habenula vs thalamus vs non-neuronal markers

celltype_order = c('MHb.1','MHb.1.2','MHb.2', 'MHb.3',
'LHb.2.7','LHb.1','LHb.1.3','LHb.1.3.4','LHb.4',
'Inhib.Thal','Excit.Thal',
'Ependymal','Astrocyte','Microglia','Endo','Oligo','OPC')

p_bubble = get_bubble_plot_sce(multiome_sce, 
  top_markers = c(
    'OLIG2', 'PDGFRA','MOG','S100B','CLDN5', 'PECAM1', 'P2RY12','CX3CR1',
'SLC1A2', 'GFAP', 'FOXJ1', 'PIFO', #Non-neuronal markers
'GAP43','SNAP25', 'SLC17A6','SLC17A7', #Neuronal markers
'LYPD6B','ADARB2','SOX14','RORB', #THalamus markers
'GAD1','GAD2','SLC32A1', #Inhibitory markers
'LYNX1','CHRM3','GABRA1','PCDH10','HTR2C', #Lateral habenula markers
'GPR151','POU4F1', #Habenula markers
'CHRNB4', 'CHAT','SLC5A7','SLC18A3','TAC1','TAC3'), #medial habenula markers
 sample_name = "Multiome Habenula: Broad marker annotations", group_col = "mid_cluster_adj", group_order = celltype_order)

p_bubble[[1]]
p_bubble[[2]] 


ggsave(p_bubble[[1]], filename = paste0(plot_path, '/all_celltypes_and_markers_meanExp.pdf'), device = 'pdf', 
width = 10, height = 6)
ggsave(p_bubble[[2]], filename = paste0(plot_path, '/all_celltypes_and_markers_zscore.pdf'), device = 'pdf', 
width = 10, height = 6)

#Some violin plots across cell-types of specific genes of interest

multiome_sce$class_label = 'Non-neurons'
multiome_sce$class_label[multiome_sce$mid_cluster_adj %in% c('MHb.1','MHb.1.2','MHb.2', 'MHb.3')] = 'MHb'
multiome_sce$class_label[multiome_sce$mid_cluster_adj %in% c('LHb.2.7','LHb.1','LHb.1.3','LHb.1.3.4','LHb.4')] = 'LHb'
multiome_sce$class_label[multiome_sce$mid_cluster_adj %in% c('Inhib.Thal','Excit.Thal')] = 'Thalamus'


target_gene = 'POU4F1'

gene_exp = assay(multiome_sce, 'scaledata')[target_gene, ]
cell_type_labels = multiome_sce$class_label
donors = multiome_sce$orig.ident
violin_df = data.frame(CPM = gene_exp, Celltype = cell_type_labels, donor = donors)

ggplot(violin_df, aes(x = Celltype, y = CPM, fill = Celltype)) +
  geom_violin(scale = 'width') +
  scale_fill_manual(values = my_colors_class) +
  labs(title = paste("Expression of", target_gene), x = "Cell Type", y = "z-score EXP") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~ donor) 






#############################
#
#
#UMAPs of the zebrafish, mouse, and yalcinbas datasets
#Will be small, and will need to try out different colors, having them all be the same would be odd
#
#
###############################


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


hu_mu_zf_ortholog_df = data.table::fread(path_to_orthologs)
table(hu_mu_zf_ortholog_df$`Mouse homology type`)

#Get the 1to1 orthologs for the mouse
hu_mu_zf_ortholog_df <- hu_mu_zf_ortholog_df %>%
  filter(
    `Mouse homology type` == 'ortholog_one2one'
  ) %>%
  filter(!duplicated(`Gene name`))
dim(hu_mu_zf_ortholog_df)

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
hashikawa_sce_sub = sce_mouse_sub$all
hashikawa_sce_neuron = sce_mouse_sub$neuron #This one has the annotated habenula clusters from the authors
rm(sce_mouse_sub)
#Gene symbols as the rownames
rownames(hashikawa_sce_sub) = rowData(hashikawa_sce_sub)$Symbol
rownames(hashikawa_sce_neuron) = rowData(hashikawa_sce_neuron)$Symbol
#Add CPM
assay(hashikawa_sce_sub, "cpm") = MetaMarkers::convert_to_cpm(assay(hashikawa_sce_sub, "counts"))
assay(hashikawa_sce_neuron, "cpm") = MetaMarkers::convert_to_cpm(assay(hashikawa_sce_neuron, "counts"))

############
#Check out OPRM1 expression using the original author annotations
############

#Use the author_subHab_celltype for this dataset

celltype_order = c('Not habenula', 'Mhb_Dorsal', 'Mhb_Superior', 'Mhb_Lateral', 'Mhb_Ventral 2/3',
'Mhb_Ventrolateral','Lhb_Hox','Lhb_Lateral','Lhb_Marginal','Lhb_Oval/Medial')

p_bubble = get_bubble_plot_sce(all_mouse_sce, 
  top_markers = c(
  'Lynx1','Chrm3','Gabra1','Pcdh10','Htr2c', #Lateral habenula markers
  'Gpr151','Pou4f1','Tac2', #Habenula markers
  'Chrnb4', 'Chat','Slc5a7','Slc18a3','Tac1'), #medial habenula markers,
 sample_name = "Wallace mouse Habenula", group_col = "author_subHab_celltype", group_order = celltype_order, exp_assay = 'cpm')
p_bubble

#Use the celltype metadata for this dataset
hashikawa_sce_neuron$celltype


celltype_order = c('MHb1','MHb5','MHb4','MHb2','MHb3','MHb6','LHb6','LHb5','LHb4','LHb3','LHb2','LHb1')

p_bubble = get_bubble_plot_sce(hashikawa_sce_neuron, 
  top_markers = c(
  'Lynx1','Chrm3','Gabra1','Pcdh10','Htr2c', #Lateral habenula markers
  'Gpr151','Pou4f1','Tac2', #Habenula markers
  'Chrnb4', 'Chat','Slc5a7','Slc18a3','Tac1'), #medial habenula markers
 sample_name = "Hashikawa mouse Habenula", group_col = "celltype", group_order = celltype_order, exp_assay = 'cpm')
p_bubble



# Add the mouse metadata for wallace
current_mouse_metadata = readRDS(paste0(wallace_03_data_path, '/wallace_mouse_metaclust_celltype_annot_metadata.rds'))
colData(all_mouse_sce) = S4Vectors::DataFrame(current_mouse_metadata)

#Load embedding
umap_mat <- readRDS(paste0(wallace_03_data_path,"/wallace_mouse_umap_embeddings.rds"))
#Reorder to SCE cell order, then attach
umap_mat <- umap_mat[colnames(all_mouse_sce), , drop = FALSE]
reducedDim(all_mouse_sce, "umap") <- umap_mat


# Add the hashikawa metadata
hashikawa_metadata = readRDS(paste0(hashikawa_data_path, '/hashikawa_mouse_metaclust_celltype_annot_metadata.rds'))
neuronal_subset_metadata = readRDS(paste0(hashikawa_data_path, '/hashikawa_mouse_neuron_metaclust_celltype_annot_metadata.rds'))
colData(hashikawa_sce_sub) = S4Vectors::DataFrame(hashikawa_metadata)

#Pass in the neuronal subset annotations
index = match(rownames(neuronal_subset_metadata), colnames(hashikawa_sce_sub) )
hashikawa_sce_sub$meta_clust_celltype_annot[index] = neuronal_subset_metadata$meta_clust_celltype_annot
#These are the annotations we used to match to Wallace
hashikawa_sce_sub$author_celltype = hashikawa_sce_sub$celltype
hashikawa_sce_sub$author_celltype[index] = neuronal_subset_metadata$celltype
table(hashikawa_sce_sub$author_celltype, hashikawa_sce_sub$celltype)

#Give a generic neuron label for the cells that were not included in the neuronal subset
table(hashikawa_sce_sub$meta_clust_celltype_annot)
index = hashikawa_sce_sub$meta_clust_celltype_annot %in% c('LHb_1.1','LHb_1.2','LHb_1.3','MHb.1','MHb.2','MHb.3')
hashikawa_sce_sub$meta_clust_celltype_annot[index] = 'Neuron'

index = hashikawa_sce_sub$author_celltype %in% c('Neuron1','Neuron2','Neuron3','Neuron4','Neuron5','Neuron6', 'Neuron7','Neuron8')
hashikawa_sce_sub$author_celltype[index] = 'Neuron'

#Load embedding
umap_mat <- readRDS(paste0(hashikawa_data_path,"/hashikawa_mouse_umap_embeddings.rds"))
#Reorder to SCE cell order, then attach
umap_mat <- umap_mat[colnames(hashikawa_sce_sub), , drop = FALSE]
reducedDim(hashikawa_sce_sub, "umap") <- umap_mat


#Filter out any cells that are listed as neurons in the full dataset but not present in the annotated neuron subset I have
#cell_subset = colnames(hashikawa_sce_sub)[hashikawa_sce_sub$celltype %in% c('Neuron1','Neuron2','Neuron3','Neuron4','Neuron5','Neuron6','Neuron7', 'Neuron8')]
#cells_exclude = cell_subset[!cell_subset %in% rownames(neuron_hashikawa_metadata)]

#hashikawa_sce_sub = hashikawa_sce_sub[, !colnames(hashikawa_sce_sub) %in% cells_exclude]

#Pass on the neuron annotations
#hashikawa_sce_sub$meta_clust_celltype_annot = hashikawa_sce_sub$celltype
#index = match(rownames(neuron_hashikawa_metadata) , colnames(hashikawa_sce_sub))
#hashikawa_sce_sub$meta_clust_celltype_annot[index] = neuron_hashikawa_metadata$meta_clust_celltype_annot


#Match annotation name
zeb_sce$final_Annotations = zeb_sce$meta_clust_celltype_annot
all_mouse_sce$final_Annotations = all_mouse_sce$meta_clust_celltype_annot
hashikawa_sce_sub$final_Annotations = hashikawa_sce_sub$meta_clust_celltype_annot



library(colorspace)

# Generate n related colors around a base color
make_cluster_shades <- function(base_color, n, span = 0.35) {
  # span controls how far from the base color you go
  shifts <- seq(-span, span, length.out = n)
  vapply(
    shifts,
    \(s) if (s < 0) darken(base_color, amount = -s) else lighten(base_color, amount = s),
    character(1)
  )
}


broad_lateral_color <- my_colors_class['LHb'] 
broad_medial_color <- my_colors_class['MHb'] 
broad_nonN_color <- my_colors_class['Non-neurons'] 


#Zebrafish palette
zeb_lat_fine <- c("ventral", "ventral_2", 'inhibitory_gap43')
zeb_med_fine <- c("dorsolateral_left_subP", "dorsomedial_neuron", 'dorsomedial_right_cholinergic', 'dorsomedial_right_cholinergic_2')

zeb_lat_fine_colors <- setNames(make_cluster_shades(broad_lateral_color, length(zeb_lat_fine)), zeb_lat_fine)
zeb_med_fine_colors <- setNames(make_cluster_shades(broad_medial_color, length(zeb_med_fine)), zeb_med_fine)

zeb_palette_vec <- c(broad_nonN_color, 'outliers' = 'grey50', zeb_lat_fine_colors, zeb_med_fine_colors)

#yalcinbas palette
yalcinbas_lat_fine <- c("LHb.1", "LHb.2", 'LHb.3', 'LHb.4', 'LHb.5', 'LHb.6', 'LHb.7')
yalcinbas_med_fine <- c('MHb.1','MHb.2','MHb.3')

yalcinbas_lat_fine_colors <- setNames(make_cluster_shades(broad_lateral_color, length(yalcinbas_lat_fine)), yalcinbas_lat_fine)
yalcinbas_med_fine_colors <- setNames(make_cluster_shades(broad_medial_color, length(yalcinbas_med_fine)), yalcinbas_med_fine)

yalcinbas_palette_vec <- c(broad_nonN_color, 'outliers' = 'grey50', 
yalcinbas_lat_fine_colors, yalcinbas_med_fine_colors,
'Thalamus' = "#670289ff")

#wallace palette
wallace_lat_fine <- c("LHb_1", "LHb_2")
wallace_med_fine <- c('MHb_cholinergic','MHb_subP','MHb_subP_cholinergic')

wallace_lat_fine_colors <- setNames(make_cluster_shades(broad_lateral_color, length(wallace_lat_fine)), wallace_lat_fine)
wallace_med_fine_colors <- setNames(make_cluster_shades(broad_medial_color, length(wallace_med_fine)), wallace_med_fine)

wallace_palette_vec <- c(broad_nonN_color, 'outliers' = 'grey50', wallace_lat_fine_colors, wallace_med_fine_colors)


#hashikawa palette
hashikawa_lat_fine <- c("LHb_1_1", "LHb_1_2", 'LHb_1_3', 'LHb_1_4', 'LHb_1_5', 'LHb_2', 'Neuron')
hashikawa_med_fine <- c('MHb_cholinergic_1','MHb_cholinergic_2','MHb_cholinergic_3', 'MHb_cholinergic_4','MHb_subP', 'MHb_subP_cholinergic')

hashikawa_lat_fine_colors <- setNames(make_cluster_shades(broad_lateral_color, length(hashikawa_lat_fine)), hashikawa_lat_fine)
hashikawa_med_fine_colors <- setNames(make_cluster_shades(broad_medial_color, length(hashikawa_med_fine)), hashikawa_med_fine)

hashikawa_palette_vec <- c(broad_nonN_color, 'outliers' = 'grey50', hashikawa_lat_fine_colors, hashikawa_med_fine_colors)



#Group the non-neurons
zeb_sce$plot_annot = zeb_sce$final_Annotations
zeb_sce$plot_annot[zeb_sce$final_Annotations %in% c('non_neuronal')] = 'Non-neurons'
zeb_sce$plot_annot[zeb_sce$final_Annotations == 'dorsomedial_right_cholinergic_GAT1'] = 'dorsomedial_right_cholinergic_2'
zeb_sce$plot_annot[zeb_sce$final_Annotations == 'dorsolateral_left_subP_BDNF'] = 'dorsolateral_left_subP'
zeb_sce$plot_annot[zeb_sce$final_Annotations == 'ventral_immediate_early'] = 'ventral_2'

yalcinbas_sce$plot_annot = yalcinbas_sce$final_Annotations
yalcinbas_sce$plot_annot[yalcinbas_sce$final_Annotations %in% c('Astrocyte', 'Endo', 'Microglia','Oligo','OPC')] = 'Non-neurons'
yalcinbas_sce$plot_annot[yalcinbas_sce$final_Annotations %in% c('Excit.Thal', 'Inhib.Thal')] = 'Thalamus'

all_mouse_sce$plot_annot = all_mouse_sce$final_Annotations
all_mouse_sce$plot_annot[all_mouse_sce$final_Annotations %in% c('Astrocytes', 'Differentiating Oligodendrocytes', 'Endothelial', 'Fibroblasts', 
'Macrophages', 'Microglia', 'Oligodendrocytes', 'Pericytes', 'Polydendrocytes')] = 'Non-neurons'

hashikawa_sce_sub$plot_annot = hashikawa_sce_sub$final_Annotations
hashikawa_sce_sub$plot_annot[hashikawa_sce_sub$final_Annotations %in% c('Astrocyte1', 'Astrocyte2','Endothelial','Epen','Microglia',
'Mural', 'Oligo1','Oligo2','Oligo3','OPC1', 'OPC2', 'OPC3')] = 'Non-neurons'


#Metacluster annots
umap_wallace_meta <- plotReducedDim(all_mouse_sce, 
                dimred = "umap",
                colour_by = "meta_cluster", 
              point_size = .5) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(
    title = "Wallace et. al UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Meta-cluster"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

umap_wallace_meta

umap_wallace <- plotReducedDim(all_mouse_sce, 
                dimred = "umap",
                colour_by = "plot_annot", 
              point_size = .5) +
  scale_color_manual(values = wallace_palette_vec) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(
    title = "Wallace et. al UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Cell type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

umap_wallace

# Add the bubble plot for how the metaclusters were annotated 
#And now the bubble plot with the Wallace marker panel

custom_markers = c('Cldn5','Abcc9','Pdgfrb','Col3a1','Cx3cr1','Mrc1','Slc6a11','Gpr17',
'Mbp','Mog','Olig1','Pdgfra','Slc17a6','Snap25','Tac2','Slc17a7','Pou4f1','Gpr151','Tac1',
'Chat', 'Slc18a3','Slc5a7', 'Gap43', 'Gad1','Gad2','Slc32a1')

cell_type_order = c('meta_cluster4', 'meta_cluster14', 'meta_cluster8','meta_cluster10',
'meta_cluster1','meta_cluster2', 'meta_cluster9', 'meta_cluster11', 'meta_cluster6',
'meta_cluster5','meta_cluster13','meta_cluster7', 'meta_cluster12', 'meta_cluster3')

p_bubble_paper_markers = get_bubble_plot_sce(all_mouse_sce , custom_markers, 
  'Wallace 2019: mouse', group_col = 'meta_cluster', group_order = cell_type_order)
p_bubble_paper_markers

ggsave(p_bubble_paper_markers[[2]], filename = paste0(plot_path, '/wallace_markers_for_meta_annot_bubble.pdf'), device = 'pdf', 
width = 8, height = 4)



umap_hashikawa_author <- plotReducedDim(hashikawa_sce_sub, 
                dimred = "umap",
                colour_by = "author_celltype", 
              point_size = .5) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(
    title = "Hashikawa et. al UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "AuthorCell type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

umap_hashikawa_author


umap_hashikawa <- plotReducedDim(hashikawa_sce_sub, 
                dimred = "umap",
                colour_by = "plot_annot", 
              point_size = .5) +
  scale_color_manual(values = hashikawa_palette_vec) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(
    title = "Hashikawa et. al UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Cell type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

umap_hashikawa

custom_markers = c('Cldn5','Foxj1','Abcc9','Pdgfrb','Col3a1','Cx3cr1','Mrc1','Slc6a11','Gpr17',
'Mbp','Mog','Olig1','Pdgfra','Slc17a6','Snap25','Tac2','Slc17a7','Pou4f1','Gpr151','Tac1',
'Chat', 'Slc18a3','Slc5a7', 'Gap43', 'Gad1','Gad2','Slc32a1')

cell_type_order = c('LHb3','LHb1','LHb2','LHb4','LHb5','LHb6','MHb6','MHb4', 'MHb3', 'MHb2','MHb1',
'MHb5','Neuron','Oligo1','Oligo2','Oligo3', 'OPC1','OPC2','OPC3', 'Astrocyte1', 'Astrocyte2','Microglia',
'Mural','Epen', 'Endothelial')

p_bubble_paper_markers = get_bubble_plot_sce(hashikawa_sce_sub , custom_markers, 
  'Hashikawa 2020: mouse', group_col = 'author_celltype', group_order = cell_type_order)
p_bubble_paper_markers

ggsave(p_bubble_paper_markers[[2]], filename = paste0(plot_path, '/hashikawa_markers_for_meta_annot_bubble.pdf'), device = 'pdf', 
width = 8, height = 4)






umap_yalcinbas <- plotReducedDim(yalcinbas_sce, 
                dimred = "UMAP",
                colour_by = "plot_annot", 
              point_size = .5) +
  scale_color_manual(values = yalcinbas_palette_vec) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(
    title = "Yalcinbas et. al UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Cell type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

umap_yalcinbas





umap_zebrafish <- plotReducedDim(zeb_sce, 
                dimred = "UMAP",
                colour_by = "plot_annot",
              point_size = .5) +
  scale_color_manual(values = zeb_palette_vec) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(
    title = "Pendey et. al UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Cell type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

umap_zebrafish_meta <- plotReducedDim(zeb_sce, 
                dimred = "UMAP",
                colour_by = "meta_cluster",
              point_size = .5) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(
    title = "Pendey et. al UMAP",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Cell type"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )


#Zebrafish regional markers
custom_markers = c('CLDN5','MRC1','SLC6A11',
'MBP','SLC17A6','SNAP25','TAC3','SLC17A7','POU4F1','GPR151','TAC1',
'SLC5A7', 'GAP43', 'GAD1','GAD2','SLC32A1')

cell_type_order = c('meta_cluster7','meta_cluster2','meta_cluster1','meta_cluster5','meta_cluster4','meta_cluster3','meta_cluster6','meta_cluster8')

p_bubble_paper_markers = get_bubble_plot_sce(zeb_sce , custom_markers, 
  'Pandey 2018: zebrafish', group_col = 'meta_cluster', group_order = cell_type_order)
p_bubble_paper_markers

ggsave(p_bubble_paper_markers[[2]], filename = paste0(plot_path, '/zebrafish_markers_for_meta_annot_bubble.pdf'), device = 'pdf', 
width = 8, height = 4)




umap_zebrafish
umap_zebrafish_meta
umap_yalcinbas
umap_wallace
umap_hashikawa

ggsave(umap_zebrafish, filename = 'zebrafish_grouped_small_umap.pdf', path = plot_path, device = 'pdf', height = 3, width = 6)
ggsave(umap_zebrafish_meta, filename = 'zebrafish_meta_cluster_umap.pdf', path = plot_path, device = 'pdf', height = 3, width = 6)

ggsave(umap_yalcinbas, filename = 'yalcinbas_grouped_small_umap.pdf', path = plot_path, device = 'pdf', height = 3, width = 6)

ggsave(umap_wallace, filename = 'wallace_grouped_small_umap.pdf', path = plot_path, device = 'pdf', height = 3, width = 6)
ggsave(umap_wallace_meta, filename = 'wallace_meta_cluster_umap.pdf', path = plot_path, device = 'pdf', height = 3, width = 6)

ggsave(umap_hashikawa, filename = 'hashikawa_grouped_small_umap.pdf', path = plot_path, device = 'pdf', height = 3, width = 6)
ggsave(umap_hashikawa_author, filename = 'hashikawa_author_annot_umap.pdf', path = plot_path, device = 'pdf', height = 3, width = 6)



#Violin plot of OPRM1 expression across the clusters

# Extract UMAP coordinates
umap_coords <- reducedDim(all_mouse_sce, "umap")

# Extract Oprm1 expression (normalized counts)
oprm1_expr <- assay(all_mouse_sce, "cpm")["Oprm1", ]

# Create dataframe with UMAP coordinates and Oprm1 expression
oprm1_umap_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  Oprm1 = oprm1_expr,
  cell_id = colnames(all_mouse_sce)
)

# Z-scale the Oprm1 expression
oprm1_umap_df$Oprm1_zscore <- scale(oprm1_umap_df$Oprm1)[, 1]
oprm1_umap_df$alpha_vec = ifelse(oprm1_umap_df$Oprm1_zscore > 0, 1, 0.5)
umap_oprm1_wallace = ggplot(oprm1_umap_df, aes(x = UMAP1, y = UMAP2, color = Oprm1_zscore, alpha = alpha_vec)) +
  geom_point(size = 0.5) +
  scale_color_gradientn(colors = c("grey", "#f48d8d", "#ff0000"), values = scales::rescale(c(0, 0.5, 1))) +
  labs(title = "Wallace mouse", x = "UMAP 1", y = "UMAP 2", color = "Oprm1 z-score") +
  guides(alpha = "none") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )


# Extract UMAP coordinates
umap_coords <- reducedDim(hashikawa_sce_sub, "umap")
oprm1_expr <- assay(hashikawa_sce_sub, "cpm")["Oprm1", ]
oprm1_umap_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  Oprm1 = oprm1_expr,
  cell_id = colnames(hashikawa_sce_sub)
)

oprm1_umap_df$Oprm1_zscore <- scale(oprm1_umap_df$Oprm1)[, 1]
oprm1_umap_df$alpha_vec = ifelse(oprm1_umap_df$Oprm1_zscore > 0, 1, 0.5)
umap_oprm1_hashikawa = ggplot(oprm1_umap_df, aes(x = UMAP1, y = UMAP2, color = Oprm1_zscore, alpha = alpha_vec)) +
  geom_point(size = 0.5) +
  scale_color_gradientn(colors = c("grey", "#f48d8d", "#ff0000"), values = scales::rescale(c(0, 0.5, 1))) +
  labs(title = "Hashikawa mouse", x = "UMAP 1", y = "UMAP 2", color = "Oprm1 z-score") +
  guides(alpha = "none") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )



# Extract UMAP coordinates
umap_coords <- reducedDim(zeb_sce, "UMAP")
oprm1_expr <- assay(zeb_sce, "cpm")["OPRM1", ]
oprm1_umap_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  Oprm1 = oprm1_expr,
  cell_id = colnames(zeb_sce)
)

oprm1_umap_df$Oprm1_zscore <- scale(oprm1_umap_df$Oprm1)[, 1]
oprm1_umap_df$alpha_vec = ifelse(oprm1_umap_df$Oprm1_zscore > 0, 1, 0.5)
umap_oprm1_zeb = ggplot(oprm1_umap_df, aes(x = UMAP1, y = UMAP2, color = Oprm1_zscore, alpha = alpha_vec)) +
  geom_point(size = 0.5) +
  scale_color_gradientn(colors = c("grey", "#f48d8d", "#ff0000"), values = scales::rescale(c(0, 0.5, 1))) +
  labs(title = "Pendey Zebrafish", x = "UMAP 1", y = "UMAP 2", color = "Oprm1 z-score") +
  guides(alpha = "none") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

# Extract UMAP coordinates
umap_coords <- reducedDim(yalcinbas_sce, "UMAP")
oprm1_expr <- assay(yalcinbas_sce, "cpm")["OPRM1", ]
oprm1_umap_df <- data.frame(
  UMAP1 = umap_coords[, 1],
  UMAP2 = umap_coords[, 2],
  Oprm1 = oprm1_expr,
  cell_id = colnames(yalcinbas_sce)
)

oprm1_umap_df$Oprm1_zscore <- scale(oprm1_umap_df$Oprm1)[, 1]
oprm1_umap_df$alpha_vec = ifelse(oprm1_umap_df$Oprm1_zscore > 0, 1, 0.5)
umap_oprm1_yalcinbas = ggplot(oprm1_umap_df, aes(x = UMAP1, y = UMAP2, color = Oprm1_zscore, alpha = alpha_vec)) +
  geom_point(size = 0.5) +
  scale_color_gradientn(colors = c("grey", "#f48d8d", "#ff0000"), values = scales::rescale(c(0, 0.5, 1))) +
  labs(title = "Yalcinbas Human", x = "UMAP 1", y = "UMAP 2", color = "Oprm1 z-score") +
  guides(alpha = "none") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )


ggsave(umap_oprm1_zeb, filename = 'oprm1_zebrafish_grouped_small_umap.pdf', path = plot_path, device = 'pdf', height = 3, width = 6)
ggsave(umap_oprm1_yalcinbas, filename = 'oprm1_yalcinbas_grouped_small_umap.pdf', path = plot_path, device = 'pdf', height = 3, width = 6)
ggsave(umap_oprm1_wallace, filename = 'oprm1_wallace_grouped_small_umap.pdf', path = plot_path, device = 'pdf', height = 3, width = 6)
ggsave(umap_oprm1_hashikawa, filename = 'oprm1_hashikawa_grouped_small_umap.pdf', path = plot_path, device = 'pdf', height = 3, width = 6)







###############################
#
#
#Select marker gene panels for the conserved cell-types across the species
#
#
###############################

#GPR151, POU4F1, CHAT, SLC5A7, SLC18A3, TAC1, TAC3/2 GAP43, GAD1, GAD2, SLC32A1

#Wallace data

all_mouse_sce$grouped_annot = all_mouse_sce$final_Annotations
all_mouse_sce$grouped_annot[all_mouse_sce$final_Annotations %in% c('MHb_subP')] = 'MHb_A'
all_mouse_sce$grouped_annot[all_mouse_sce$final_Annotations %in% c('MHb_cholinergic')] = 'MHb_B'
all_mouse_sce$grouped_annot[all_mouse_sce$final_Annotations %in% c('LHb_2')] = 'LHb_A'
all_mouse_sce$grouped_annot[all_mouse_sce$final_Annotations %in% c('LHb_1')] = 'LHb_C'

celltype_order = c('LHb_C','LHb_B','LHb_A','MHb_B','MHb_subP_cholinergic','MHb_A')
cell_index = all_mouse_sce$grouped_annot %in% celltype_order

p_bubble = get_bubble_plot_sce(all_mouse_sce[ , cell_index], 
  top_markers = c('Tac1','Tac2',
  'Sstr2','C1ql1',
  'Chat','Slc5a7', 'Slc18a3',
  'Oprm1', 'Pcdh10', 'Htr2c',
  'Col25a1',
  'Grik4', 
  'Gap43','Rph3a','Chrm3', 
  'Gad1','Gad2','Slc32a1',
  'Gpr151','Pou4f1'),
 sample_name = "Wallace Mouse", group_col = "grouped_annot", group_order = celltype_order)

#p_bubble[[1]]
wallace_bubble = p_bubble[[2]] 
wallace_bubble


#Hashikawa data
hashikawa_sce_sub$grouped_annot = hashikawa_sce_sub$final_Annotations
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$final_Annotations %in% c('MHb_subP')] = 'MHb_A'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$final_Annotations %in% c(paste0('MHb_cholinergic_', 1:4))] = 'MHb_B'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$final_Annotations %in% c('LHb_2', 'LHb_1_5')] = 'LHb_A'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$final_Annotations %in% c('LHb_1_2', 'LHb_1_4')] = 'LHb_B'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$final_Annotations %in% c('LHb_1_1', 'LHb_1_3')] = 'LHb_C'

celltype_order = c('LHb_C','LHb_B','LHb_A','MHb_B','MHb_subP_cholinergic','MHb_A')
cell_index = hashikawa_sce_sub$grouped_annot %in% celltype_order

p_bubble = get_bubble_plot_sce(hashikawa_sce_sub[ , cell_index], 
  top_markers = c('Tac1','Tac2',
  'Sstr2','C1ql1',
  'Chat','Slc5a7', 'Slc18a3',
  'Oprm1', 'Pcdh10', 'Htr2c',
  'Col25a1',
  'Grik4', 
  'Gap43','Rph3a','Chrm3', 
  'Gad1','Gad2','Slc32a1',
  'Gpr151','Pou4f1'),
 sample_name = "Hashikawa Mouse", group_col = "grouped_annot", group_order = celltype_order)

#p_bubble[[1]]
hashikawa_bubble = p_bubble[[2]] 
hashikawa_bubble 


#Pandey data

zeb_sce$grouped_annot = zeb_sce$plot_annot
zeb_sce$grouped_annot[zeb_sce$plot_annot %in% c('dorsolateral_left_subP')] = 'MHb_A'
zeb_sce$grouped_annot[zeb_sce$plot_annot %in% c('dorsomedial_right_cholinergic','dorsomedial_right_cholinergic_2')] = 'MHb_B'
zeb_sce$grouped_annot[zeb_sce$plot_annot %in% c('ventral', 'ventral_2')] = 'LHb_B'
zeb_sce$grouped_annot[zeb_sce$plot_annot %in% c('inhibitory_gap43')] = 'LHb_C'


cell_index = !zeb_sce$grouped_annot %in% c('outliers', 'Non-neurons')

celltype_order = c('LHb_C','LHb_B','LHb_A','MHb_B','dorsomedial_neuron','MHb_A')

p_bubble = get_bubble_plot_sce(zeb_sce[ , cell_index], 
  top_markers = c('TAC1','TAC3',
  'SLC5A7', 'SLC18A3',
  'OPRM1','PCDH10','HTR2C',  
  'GRIK4', 
  'GAP43','CHRM3', 
  'GAD1','GAD2','SLC32A1',
  'GPR151','POU4F1'),
 sample_name = "Pandey Zebrafish", group_col = "grouped_annot", group_order = celltype_order)

#p_bubble[[1]]
pandey_bubble = p_bubble[[2]] 

pandey_bubble


# #Yalcinbas data
# yalcinbas_sce$grouped_annot = yalcinbas_sce$final_Annotations
# yalcinbas_sce$grouped_annot[yalcinbas_sce$final_Annotations %in% c('MHb.1')] = 'MHb.1'
# yalcinbas_sce$grouped_annot[yalcinbas_sce$final_Annotations %in% c('LHb.6')] = 'MHb.1.2'
# yalcinbas_sce$grouped_annot[yalcinbas_sce$final_Annotations %in% c('MHb.2')] = 'MHb.2'
# yalcinbas_sce$grouped_annot[yalcinbas_sce$final_Annotations %in% c('LHb.2', 'LHb.7')] = 'LHb.2.7'
# yalcinbas_sce$grouped_annot[yalcinbas_sce$final_Annotations %in% c('LHb.5', 'LHb.1')] = 'LHb.1.3.4'
# yalcinbas_sce$grouped_annot[yalcinbas_sce$final_Annotations %in% c('LHb.3', 'LHb.4')] = 'LHb.4'

# celltype_order = c('LHb.4','LHb.1.3.4','LHb.2.7','MHb.2','MHb.1.2','MHb.1')
# cell_index = yalcinbas_sce$grouped_annot %in% celltype_order

# p_bubble = get_bubble_plot_sce(yalcinbas_sce[ , cell_index], 
#   top_markers = c('TAC1','TAC3',
#   'SSTR2','C1QL1',
#   'CHAT','SLC5A7', 'SLC18A3',
#   'OPRM1','PCDH10','HTR2C', 
#   'COL25A1', 
#   'GRIK4', 
#   'GAP43','RPH3A','CHRM3', 
#   'GAD1','GAD2','SLC32A1',
#   'GPR151','POU4F1'),
#  sample_name = "Yalcinbas Human", group_col = "grouped_annot", group_order = celltype_order)

# #p_bubble[[1]]
# yalcinbas_bubble = p_bubble[[2]] 
# yalcinbas_bubble

# multiome_sce$final_Annotations = as.character(multiome_sce$mid_cluster_adj)
# multiome_sce$final_Annotations[ multiome_sce$refined_mid_cluster == 'Inhib_LHb_4.1'] = 'Inhib_LHb_4.1'
# multiome_sce$final_Annotations[ multiome_sce$refined_mid_cluster == 'Inhib_LHb_4.2'] = 'Inhib_LHb_4.2'

# multiome_sce$grouped_annot = multiome_sce$final_Annotations
# multiome_sce$grouped_annot[multiome_sce$final_Annotations %in% c('LHb.1', 'LHb.1.3', 'LHb.1.3.4')] = 'LHb.1.3.4'

# celltype_order = c('Inhib_LHb_4.2', 'Inhib_LHb_4.1','LHb.4','LHb.1.3.4','LHb.1.3','LHb.1','LHb.2.7','MHb.2','MHb.1.2','MHb.1')
# cell_index = multiome_sce$grouped_annot %in% celltype_order

# p_bubble = get_bubble_plot_sce(multiome_sce[ , cell_index], 
#   top_markers = c('TAC1','TAC3',
#   'SSTR2','C1QL1',
#   'CHAT','SLC5A7', 'SLC18A3',
#   'OPRM1','PCDH10','HTR2C', 
#   'COL25A1', 
#   'GRIK4', 
#   'GAP43','RPH3A','CHRM3', 
#   'GAD1','GAD2','SLC32A1',
#   'GPR151','POU4F1'),
#  sample_name = "Multiome Human", group_col = "grouped_annot", group_order = celltype_order)

# #p_bubble[[1]]
# multiome_bubble = p_bubble[[2]] 
# multiome_bubble 


wallace_bubble
hashikawa_bubble
pandey_bubble
#yalcinbas_bubble
#multiome_bubble

ggsave(wallace_bubble, filename = 'wallace_grouped_bubble.pdf', path = plot_path, device = 'pdf', height = 4, width = 6)
ggsave(hashikawa_bubble, filename = 'hashikawa_grouped_bubble.pdf', path = plot_path, device = 'pdf', height = 4, width = 6)
ggsave(pandey_bubble, filename = 'pandey_grouped_bubble.pdf', path = plot_path, device = 'pdf', height = 4, width = 8)
#ggsave(yalcinbas_bubble, filename = 'yalcinbas_grouped_bubble.pdf', path = plot_path, device = 'pdf', height = 4, width = 6)
#ggsave(multiome_bubble, filename = 'multiome_grouped_bubble.pdf', path = plot_path, device = 'pdf', height = 4, width = 6)



#Pull together non-human metadata for supplemental tables
all_mouse_sce$cell_id = colnames(all_mouse_sce)
wallace_meta = colData(all_mouse_sce) %>% as.data.frame() %>% dplyr::select(c('cell_id', 'author_celltype', 'author_subHab_celltype', 
'meta_cluster', 'meta_clust_celltype_annot','grouped_annot'))

colnames(wallace_meta) = c('cell_id','author_celltype','author_subHab_celltype','meta_cluster','Annotated_meta_cluster','consensus_annotation')


hashikawa_sce_sub$cell_id = colnames(hashikawa_sce_sub)
hashikawa_meta = colData(hashikawa_sce_sub) %>% as.data.frame() %>% dplyr::select(c('cell_id', 'author_celltype', 
'meta_clust_celltype_annot','grouped_annot'))

colnames(hashikawa_meta) = c('cell_id','author_celltype','Annotated_meta_cluster','consensus_annotation')

zeb_sce$cell_id = colnames(zeb_sce)
zeb_meta = colData(zeb_sce) %>% as.data.frame() %>% dplyr::select(c('cell_id', 'meta_cluster', 'plot_annot', 'grouped_annot'))

colnames(zeb_meta) = c('cell_id','meta_cluster','Annotated_meta_cluster','consensus_annotation')


write.csv(wallace_meta, file = paste0(new_data_path, '/wallace_meta_for_supplemental_table.csv'), row.names = F)
write.csv(hashikawa_meta, file = paste0(new_data_path, '/hashikawa_meta_for_supplemental_table.csv'), row.names = F)
write.csv(zeb_meta, file = paste0(new_data_path, '/zeb_meta_for_supplemental_table.csv'), row.names = F)

