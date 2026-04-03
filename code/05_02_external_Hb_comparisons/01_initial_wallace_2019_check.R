#Script to check the initial Habenula data from Wallace 2019

library(Seurat)

hab_batch1_data = readRDS('processed-data/05_02_external_Hb_comparisons/Wallace_etal_2019_habenula_scseq/hab_batch1.rds')
hab_batch1_data


cell_metadata = hab_batch1_data@meta.data

#Provided metadata does not contain cell type annotations, sad :(
View(cell_metadata)
dim(cell_metadata)


hab_csv_data_df = data.table::fread('processed-data/05_02_external_Hb_comparisons/Wallace_etal_2019_habenula_scseq/habData.csv')
View(hab_csv_data_df)
dim(hab_csv_data_df)



