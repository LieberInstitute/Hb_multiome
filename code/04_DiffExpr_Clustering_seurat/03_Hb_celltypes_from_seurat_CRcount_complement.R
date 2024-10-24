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
## Note: Need 80 GB mem for 100K cells
##
## Authors. CSC 
## Date. Oct 23rd 2024
########################################################################

## load libraries
library(Seurat)
library(tidyverse)
library(dplyr)
library(here)

here::here()

# Check/create directories

inputDir <- here("processed-data", "03_pseudobulking", "cellranger_count")
outputDir <- here("processed-data", "03_pseudobulking", "cellranger_count")
cross_barcodes_Dir <- here("processed-data", "100_cell_match_cellranger_gex_atac", "reanalize_files")

## Read the CSV files with the cross-barcodes to be removed from each sample

barcodes_paths <- list.files(cross_barcodes_Dir, pattern = "[0-9]+C[_-]Hb", full.names = TRUE)
basename(barcodes_paths)
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
count_mtx_type <- 'norm_counts' 
Seurat_reduction <- 'Harmony' 
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

## Load Seurat integrated with cluster information
SeuratOBJ <- get_seurat(here(inputDir, paste0(Seurat_base_name, ".rds")))

## Build df with corresponding samples IDs and CSV files with cross-barcodes to remove

df_sampleIDs_CSV <- data.frame(path=barcodes_paths, cell_prefix=paste0(parse_number(basename(barcodes_paths)), "C_"))
## Get corresponding orig.ident and cell_prefix
sample_IDs <- unique(SeuratOBJ@meta.data$orig.ident)
# [1] "10C_Hb_KDM" "11C_Hb_KDM" "12C_Hb_KDM" "3C_Hb_KDM"  "4C_Hb_KDM" 
# [6] "5C_Hb_KDM"  "6C_Hb_KDM"  "7C_Hb_KDM"  "8C_Hb_KDM"  "9C_Hb_KDM" 
cell_prefix <- paste0(unlist(readr::parse_number(sample_IDs)), "C_")
## Build df with corresponding IDs from Seurat object
df_sampleIDs_Seurat <- data.frame(orig.ident=sample_IDs, cell_prefix=cell_prefix)
## Merge table to re-format cells to be removed 
df_sampleIDs_Seurat <- left_join(df_sampleIDs_Seurat, df_sampleIDs_CSV, by="cell_prefix")

## Get some stats for validation before remove cells
print(table(SeuratOBJ$orig.ident))
# 10C_Hb_KDM 11C_Hb_KDM 12C_Hb_KDM  3C_Hb_KDM  4C_Hb_KDM  5C_Hb_KDM  6C_Hb_KDM 
# 12303      14591       5616       6275       7770       3775       8852 
# 7C_Hb_KDM  8C_Hb_KDM  9C_Hb_KDM 
# 10003       6161       8831 
table(Idents(SeuratOBJ))
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 7723 6460 5000 4775 4662 3938 3791 3614 3552 3551 3539 3412 3154 2795 2509 2475 
# 16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31 
# 2419 2074 1998 1989 1928 1771 1343 1308 1241  988  914  498  265  182  166  143 


message("Starting process to remove cross-barcodes from Cell Ranger-count Seurat: `", Seurat_base_name, "`")

#head(SeuratOBJ@active.ident)
#SeuratOBJ@meta.data[SeuratOBJ@meta.data$orig.ident==df_sampleIDs_Seurat$orig.ident[4],]

for (r in 1:nrow(df_sampleIDs_Seurat)) {
  message("\nSample: ", basename(df_sampleIDs_Seurat$path[r]))
  ## read CSV with cells to remove and format them according with corresponding Seurat `cell-prefix` format
  barcodes_to_remove <- read.csv(here(df_sampleIDs_Seurat$path[r]), header = TRUE)
  barcodes_to_remove <- unlist(barcodes_to_remove)
  barcodes_to_remove <- paste0(df_sampleIDs_Seurat$cell_prefix[r], barcodes_to_remove)
  message("Cells to remove: ", length(barcodes_to_remove))
  ## remove cells from integrated seurat
  SeuratObj_subset <- subset(SeuratOBJ, cells = barcodes_to_remove, invert = TRUE)
  }

#head(barcodes_to_remove)
#SeuratObj_subset <- subset(SeuratOBJ, cells = barcodes_to_remove, invert = TRUE)

message("Resulting subset excluding cross-barcodes from Cell Ranger-count Seurat object: `", Seurat_base_name, "`")

print(table(SeuratObj_subset$orig.ident))
table(Idents(SeuratObj_subset))

## Save new Seurat-subset 
rds_name <- paste0(Seurat_base_name,'_', '_cellRanger_count_subset.rds')
rds_name <- here(processedDir, rds_name)
saveRDS(SeuratOBJ, file = rds_name)

# SeuratOBJ[ (SeuratOBJ@meta.data$orig.ident == "3C_Hb_KDM") == TRUE) ]
# Idents(SeuratOBJ)
# length(Cells(SeuratOBJ))

