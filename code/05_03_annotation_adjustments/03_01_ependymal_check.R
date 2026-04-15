#Quick check on and ependymal marker expression in the multiome data




library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
library(dplyr)
library(ggplot2)
library(scater)
library(scran)
library(qs2)
library(here)

here::here()

#Path for new data generated
new_data_path = here('processed-data', '05_03_annotation_adjustments','03_01_ependymal_check')
#Path to plot directory
plot_path = here('plots', '05_03_annotation_adjustments','03_01_ependymal_check')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)

#Source the bubble plot functions
source(here('code','05_02_external_Hb_comparisons', 'bubble_plot_functions.R'))


multiome_path = here('processed-data', '05_01_drop_doublets','01_drop_doublets_and_reDimReduce')

#Multiome human data
multiome_sce = qs_read(paste0(multiome_path, '/reprocessed_doubletRemoved_multiomeHab_SCE.qs2'))
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))


colnames(colData(multiome_sce))

#Marker plot looking at astrocyte and ependymal markers
p_bubble = get_bubble_plot_sce(multiome_sce, 
  top_markers = c('FOXJ1','TPPP3','PIFO','CD24','TMEM212','RARRES2','CLDN5','TJP1',
  'S100B', 'GFAP','AQP4','SOX9','HOPX','SLC1A2', 
  'P2RY12','TMEM119','CX3CR1'),
 sample_name = "Multiome Habenula", group_col = "cluster_ann")

p_bubble[[1]]
p_bubble[[2]]


