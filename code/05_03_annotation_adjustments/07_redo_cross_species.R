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


multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')

multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

multiome_sce


