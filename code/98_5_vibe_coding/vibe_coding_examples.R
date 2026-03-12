
#Testing out Positron's AI Assistant feature for generating different kinds of plots for single cell data

#First, orient everyone to the Positron interface
#Show how to check on memory usage



#Load libraries and data

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(here)

#Mouse habenula data
hashikawa_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data'
list.files(hashikawa_path)

#Loads as sce_mouse_sub, list of sce objects
load(paste0(hashikawa_path, '/sce_mouse_habenula.Rdata'))
hashikawa_sce = sce_mouse_sub$all
# Add cpms for later
assay(hashikawa_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(hashikawa_sce, 'counts'))


hashikawa_sce 

#Swap the rownames to gene symbols
rowData(hashikawa_sce)
rownames(hashikawa_sce) = rowData(hashikawa_sce)$Symbol

#Check out the available metadata

colnames(colData(hashikawa_sce))

table(hashikawa_sce$celltype)

table(hashikawa_sce$stim)


#Make a proportional bar plot of cell types by stimulus condition



#Adjust color palette


#Try different palettes with subjective descriptions. Try different models







#Make a bubble plot for a marker gene panel
#script to grab marker genes from: 98_external_Hb_comparisons/03_metaMarkers_wallace_2019

#Try the agent to grab the markers from that script




#Make the pubble plot



# Add a z-scored version




#Try a different model







#Compute cell-type markers

celltype_markers = MetaMarkers::compute_markers(assay(hashikawa_sce, 'cpm'), hashikawa_sce$celltype)

#View the top 10 markers by auroc




#Make a heatmap of the top 10 markers per cell type




#Highlight specific genes
# Create heatmap with Apoe highlighted




#Maybe try recommendations for dimensionality rededuction plots



#Try visualizing average marker gene expression on a UMAP


##################
# 
# 
# Any other plot suggestions?
#
#  
##################
