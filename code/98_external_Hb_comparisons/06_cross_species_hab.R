#Comparing the habenula metaclusters across zebrafish and mouse with MetaNeighbor
#Might be difficult due to the zebrafish genome duplication, a lot of informative genes have one-many matches between mouse and zebrafish
#Those one-many genes will get filtered out with the traditional approach of just keeping the one-to-ones


library(SingleCellExperiment)
library(Seurat)
library(MetaNeighbor)
library(dplyr)
library(ggplot2)
library(here)

here::here()



#Path to Wallace annotated data
wallace_02_data_path = here('processed-data', '98_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')

#Load up the full SCE object, contains the metacluster annotations
all_mouse_sce = readRDS(paste0(wallace_02_data_path, '/all_donor_sce_with_denovo_clusters.rds'))
#Convert to seurat and then add in the annotated metdata
all_mouse_seurat = as.Seurat(all_mouse_sce)



