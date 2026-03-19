
#Sub clustering the LHb.4 population
#The co-GABA-GLut population is actually very small, and this is the largest cluster we have for the habenula.
#See if we can resolve some higher resolution clusters.



library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
library(MetaNeighbor)
library(dplyr)
library(ggplot2)
library(scater)
library(scran)
library(here)

here::here()

#Path for new data generated
new_data_path = here('processed-data', '99_donor_cluster_replicability','03_LHb4_subclustering')
#Path to plot directory
plot_path = here('plots', '99_donor_cluster_replicability','03_LHb4_subclustering')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '08_spatial_registration_vs_multiome_snRNA-seq','mid')

#Multiome human data
multiome_sce = readRDS(paste0(multiome_path, '/seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v5.rds'))
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))









