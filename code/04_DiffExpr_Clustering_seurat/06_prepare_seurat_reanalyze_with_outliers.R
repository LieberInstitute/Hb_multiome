library(tidyverse)
library(dplyr)
library(data.table)
library(magrittr)
library(here)

here::here()


################## (1) Load integrated Seurat with Harmony correction


# Check/create directories
inputDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze")
inputDir_cvs <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze", "cvs_files_markers")
processedDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze")
cvsDir <- here("processed-data", "04_DiffExpr_Clustering_seurat", "cellrangerARC_reanalyze", "cvs_files_markers")

## Check directories
if (!dir.exists(processedDir)) {dir.create(processedDir)}
if (!dir.exists(cvsDir)) {dir.create(cvsDir)}

## Contains marker lists 
source(here("code", "04_DiffExpr_Clustering_seurat", "remote_DGE_marker_gene_lists.R"))       # Call functions to read paths

get_seurat <- function(name) { sobj <- readRDS(name)}


#############################           Initials        ################################

## Set count-mtx type and integration model (CCA or Harmony)

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

message("Starting cell-type identification for ", Seurat_base_name)

## Load Seurat Integrated with cluster information
SeuratOBJ <- get_seurat(here(inputDir, paste0(Seurat_base_name, ".rds")))
## verification of the integration
print(table(SeuratOBJ$orig.ident))
#SeuratOBJ@reductions


################## (2) load vector with valid barcodes by sample


## Scans arguments invoked from slurm job shell sh
sample_tmp <- commandArgs(trailingOnly = TRUE)
# For testing:
# sample_tmp <- "4S_Hb_KDM_reanalysis, 4S_Hb_KDM"
sample_data = unlist(strsplit(sample_tmp,","))
Seurat_base_name <- trimws(sample_data[[2]])

## Read barcodes after remove outliers (PASS)
csvDir_barcodes <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", "csv_files")
csvDir_barcodes <- here(csvDir_barcodes, paste0(Seurat_base_name, "_bc_PASS_isOutliers.csv"))
df_valid_barcodes_filtered <- read.csv(csvDir_barcodes)
# AAACAGCCAGAATGAC-1
# AAACAGCCAGCAAGGC-1
#head(df_valid_barcodes_filtered)
len_valid_bc <- length(df_valid_barcodes_filtered$x)
v_valid_barcodes_filtered <- df_valid_barcodes_filtered$x

message("Vector with ", len_valid_bc ," valid barcodes for sample `", Seurat_base_name,"` loaded")


################## (3) kept only cells with Outliers 


################## (4) Save Seurat object with ONLY Outliers to identify cell types with `01_Hb_celltypes_from_seurat_reanalyze.R` script



