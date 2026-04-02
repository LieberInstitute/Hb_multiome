#Consolidating the annotation changes and saving, this will be used for all downstream analyses

#Annotation changes
#C.11.MHb.1.2 changes from MHb.1.2 to MHb.2
#LHb.4 gets split into LHb.4 and inhib.LHb.4.1 and inhib.LHb.4.2
#LHb.1, LHb.1.3, and LHb.1.3.4 all get merged into LHb.1.3.4



library(SingleCellExperiment)
library(Seurat)
library(qs2)
library(dplyr)
library(ggplot2)
library(here)


here::here()

#Path for new data generated
new_data_path = here('processed-data', '99_donor_cluster_replicability','05_refined_annotations')
#Path to plot directory
plot_path = here('plots', '99_donor_cluster_replicability','05_refined_annotations')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '05_5_drop_doublets','01_drop_doublets_and_reDimReduce')

#Multiome human data, singlecellexperiment
multiome_sce = qs_read(paste0(multiome_path, '/reprocessed_doubletRemoved_multiomeHab_SCE.qs2'))

## Seurat object 
midSeurat = qs_read(paste0(multiome_path, '/reprocessed_doubletRemoved_multiomeHab_seurat.qs2'))

#This has the inhibitory LHb annotations
script_02_data_path = here('processed-data', '99_donor_cluster_replicability','02_LHb4_investigation')
multiome_seurat_integrated = readRDS(paste0(script_02_data_path, '/multiome_LHb4_LHb7_integrated_seurat.rds'))
colnames(multiome_seurat_integrated[[]])

inhib_meta = multiome_seurat_integrated[[]]

#First add the inhibitory annotations

# Add the putative Inhibitory neuron annotations to the midSeurat object
inhib_barcodes_1 = rownames(inhib_meta)[inhib_meta$refined_mid_cluster == 'Putative_Inhib_LHb_4.1']
inhib_barcodes_2 = rownames(inhib_meta)[inhib_meta$refined_mid_cluster == 'Putative_Inhib_LHb_4.2']


midSeurat$refined_mid_cluster = midSeurat$mid_cluster
midSeurat$refined_mid_cluster[rownames(midSeurat[[]]) %in% inhib_barcodes_1] = 'Putative_Inhib_LHb_4.1'
midSeurat$refined_mid_cluster[rownames(midSeurat[[]]) %in% inhib_barcodes_2] = 'Putative_Inhib_LHb_4.2'

table(midSeurat$refined_mid_cluster, midSeurat$mid_cluster)

multiome_sce$refined_mid_cluster = multiome_sce$mid_cluster
multiome_sce$refined_mid_cluster[rownames(colData(multiome_sce)) %in% inhib_barcodes_1] = 'Putative_Inhib_LHb_4.1'
multiome_sce$refined_mid_cluster[rownames(colData(multiome_sce)) %in% inhib_barcodes_2] = 'Putative_Inhib_LHb_4.2'


table(multiome_sce$refined_mid_cluster,multiome_sce$mid_cluster)



