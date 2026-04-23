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
#new_data_path = here('processed-data', '05_03_annotation_adjustments', '09_figure_plots')
#Path to plot directory
plot_path = here('plots','05_03_annotation_adjustments', '09_figure_plots')

#if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))
colnames(colData(multiome_sce))

SingleCellExperiment::altExp(multiome_sce, "ATAC") <- NULL


source(here('code','05_02_external_Hb_comparisons', 'bubble_plot_functions.R'))


#Colors
my_colors <- c(
    LHb = "#1f78b4",
    MHb = "#ad1d8c",
    Oligo = "#384a08",
    Astrocyte = "#532222", 
    OPC = "#829454",
    Microglia = "#141b02",
    Endo = "#d95f02",
    Inhib_Thal = "#9a9fe7",
    Excit_Thal = "#42467b",
    Thal = "#4d55b7",
    Ependymal = "#f5a105ff"
)

my_colors_class <- c(
    LHb = "#1f78b4",
    MHb = "#ad1d8c",
    `Non-neurons` = "#532222", 
    Thalamus = "#4d55b7"
    
)

## assign color gradients to mid resolution clusters based on Broad cell-types

# extract LHb and MHb clusters
cluster_levels <- c(unique(multiome_sce$mid_cluster))
cluster_levels
LHb_clusters <- grep("LHb", cluster_levels, value = TRUE)
MHb_clusters <- grep("MHb", cluster_levels, value = TRUE)

# Create tonal gradients for LHb and MHb
LHb_colors <- colorspace::sequential_hcl(length(LHb_clusters), h = 210, c = 80, l = c(30, 80))
MHb_colors <- colorspace::sequential_hcl(length(MHb_clusters), h = 320, c = 80, l = c(30, 80))

# Build full cluster color map
my_colors_mid <- setNames(rep("#bdbdbd", length(cluster_levels)), cluster_levels)
my_colors_mid[LHb_clusters] <- LHb_colors
my_colors_mid[MHb_clusters] <- MHb_colors

# assign base color for other types from your existing palette
for (category in c("Oligo", "Astrocyte", "OPC", "Microglia", "Endo", "Inhib.Thal", "Excit.Thal", "Thal")) {
    matched <- grep(category, cluster_levels, value = TRUE)
    my_colors_mid[matched] <- my_colors[[gsub("\\.", "_", category)]]
}

my_colors_mid["Inhib_LHb_4.1"] <- "#8B0000"  # Dark red
my_colors_mid["Inhib_LHb_4.2"] <- "#DC143C"  # Crimson red
#Adjust color for MHb3, too light
my_colors_mid["MHb.3"] <- "#56204eff" 
my_colors_mid["Inhib.Thal"] <- "#9a9fe7"
my_colors_mid["Ependymal"] <- "#f5a105ff"

names(my_colors_mid)


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

plt1 <- plotReducedDim(multiome_sce, 
                dimred = "wnn.umap",
                colour_by = "mid_cluster_adj", 
              point_size = .5) +
  scale_color_manual(values = my_colors_mid) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
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
  guides(colour = guide_legend(override.aes = list(size = 3))) +
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
  guides(colour = guide_legend(override.aes = list(size = 3))) +
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


ggsave(plt1, filename = 'wnn_umap_starting_mid_cluster.pdf', path = plot_path, device = 'pdf', 
height = 6, width = 6)
ggsave(plt2, filename = 'RNA_umap_starting_mid_cluster.pdf', path = plot_path, device = 'pdf', 
height = 6, width = 6)
ggsave(plt3, filename = 'ATAC_umap_starting_mid_cluster.pdf', path = plot_path, device = 'pdf', 
height = 6, width = 6)


#And pull together the broader marker bubble plot, here we want the habenula vs thalamus vs non-neuronal markers

celltype_order = c('MHb.1','MHb.1.2','MHb.2', 'MHb.3',
'LHb.2.7','LHb.1','LHb.1.3','LHb.1.3.4','LHb.4',
'Inhib.Thal','Excit.Thal',
'Ependymal','Astrocyte','Microglia','Endo','Oligo','OPC')

p_bubble = get_bubble_plot_sce(multiome_sce, 
  top_markers = c('OLIG2', 'PDGFRA','MOG','S100B','CLDN5', 'PECAM1', 'P2RY12','CX3CR1',
'SLC1A2', 'GFAP', 'FOXJ1', 'PIFO','GAP43','SNAP25', 'SLC17A6','SLC17A7','GAD1','GAD2','SLC32A1',
'SOX14','RORB','DRD2','LYPD6B','ADARB2',
'GPR151','POU4F1','HTR2C','CHRNB4', 'CHAT','SLC5A7','SLC18A3','TAC1','TAC3'),
 sample_name = "Multiome Habenula: Broad marker annotations", group_col = "mid_cluster_adj", group_order = celltype_order)

p_bubble[[1]]
p_bubble[[2]] 

c('OLIG2', 'PDGFRA','MOG','S100B','CLDN5', 'PECAM1', 'P2RY12','CX3CR1',
'SLC1A2', 'GFAP', 'FOXJ1', 'PIFO','GAP43','SNAP25', 'SLC17A6','SLC17A7','GAD1','GAD2','SLC32A1', 'GPR151','POU4F1')




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







