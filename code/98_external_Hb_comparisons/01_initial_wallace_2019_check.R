#Script to check the initial Habenula data from Wallace 2019

library(Seurat)

hab_batch1_data = readRDS('processed-data/98_external_Hb_comparisons/Wallace_etal_2019_habenula_scseq/hab_batch1.rds')

cell_metadata = hab_batch1_data@meta.data

#Provided metadata does not contain cell type annotations, sad :(
View(cell_metadata)





