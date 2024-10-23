########################################################################
## Obtains from the Cell Ranger-count side the complementary barcodes after the intersection between `Cell Ranger-count` and `Cell Ranger-atac`
##      for further analysis (cell-type identification and perceptual values of neu, hb and thal)
##
## INPUT:
##      Seurat objects from `Cell Ranger-count`
##      CSV file with cross-barcodes between between `Cell Ranger-count` and `Cell Ranger-atac`
## 
## OUPUT:
##      Seurat objects with the complement barcodes resulted 
##
## Authors. CSC 
## Date. Oct 23rd 2024
########################################################################

## load libraries
library(tidyverse)
library(dplyr)
library(data.table)
library(magrittr)
library(here)

here::here()

# Check/create directories
inputDir <- here("processed-data", "03_pseudobulking", "cellranger_count")
#inputDir_cvs <- here("processed-data", "03_pseudobulking", "cellranger_count", "cvs_files_markers")
#processedDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellranger_count")
#cvsDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellranger_count", "cvs_files_markers")
## Check directories
if (!dir.exists(processedDir)) {dir.create(processedDir)}
if (!dir.exists(cvsDir)) {dir.create(cvsDir)}


## Read the csv files with the barcodes to be removed from each sample

cross_barcodes_Dir <- here("processed-data", "100_cell_match_cellranger_gex_atac", "reanalize_files")
cross_bc_Dir <- here(processedDir_bc)
list.files(cross_barcodes_Dir, pattern = "[0-9]+C[_-]Hb")
# [1] "10C_Hb_KDM_matching_ARC_barcodes_1.csv"  
# [2] "11C_Hb_KDM_matching_ARC_barcodes_1.csv"  
# [3] "12C_Hb_KDM_matching_ARC_barcodes_1.csv"  
# [4] "1C-Hb-KDM-Hb_matching_ARC_barcodes_1.csv"
# [5] "2C-Hb-KDM-Hb_matching_ARC_barcodes_1.csv"
# [6] "3C_Hb_KDM_matching_ARC_barcodes_1.csv"   
# ...

## Read integrated Harmony Seurat 

get_seurat <- function(name) { sobj <- readRDS(name)}

## Define count-mtx type and integration model (CCA or Harmony)
#count_mtx_type <- 'data_counts'
count_mtx_type <- 'norm_counts' 
#Seurat_reduction <- 'CCA'
Seurat_reduction <- 'Harmony' 
## Minimum cells by cluster
minCells <- 1

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts'}

## Build Seurat object name. `subset` suffix means clusters with fewer cells than `minCells` had been filtered. 
if (Seurat_reduction=='CCA') {
  Seurat_base_name <- paste0(Seurat_base_name, '_CCA_All')
} else {
  Seurat_base_name <- paste0(Seurat_base_name, '_Harmony_All')
}
## Validate seurat exists
if (length(list.files(inputDir, pattern = Seurat_base_name)==1)) {
  message("Processing ", Seurat_base_name)
} else {
  message("Input seurat object missed!")
  stop()
}

message("Removing cross-barcodes from Cell Ranger-count Seurats `", Seurat_base_name, "`")

## Load Seurat Integrated with cluster information
SeuratOBJ <- get_seurat(here(inputDir, paste0(Seurat_base_name, ".rds")))
## verification of the integration
print(table(SeuratOBJ$orig.ident))
#SeuratOBJ@reductions



