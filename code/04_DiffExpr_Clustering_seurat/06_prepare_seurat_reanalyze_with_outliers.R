########################################################################
## Prepare a Seurat object with outliers to further identify cell-types
##  
## INPUT:
##      CSV files with barcodes that passed outliers calculated only on the GEX multiome side
## 
## OUPUT:
##      Seurat with only cells detected with Outliers 
##
## Authors. CSC 
## Date. Nov ,2024
########################################################################

library(here)

here::here()


################## (1) Load integrated Seurat with Harmony correction


# Check/create directories
inputDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze")
outputDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze_outliers") 

## Check directories
if (!dir.exists(outputDir)) {dir.create(outputDir)}

get_seurat <- function(name) { sobj <- readRDS(name)}


#############################           Initials        ################################

## Set count-mtx type and integration model (CCA or Harmony)

count_mtx_type <- 'norm_counts' 
Seurat_reduction <- 'Harmony' 
minCells <- 1

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts'}
Seurat_base_name <- paste0(Seurat_base_name, '_Harmony_All')

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
# 4S_Hb_KDM_reanalysis  5S_Hb_KDM_reanalysis  6S_Hb_KDM_reanalysis 
# 7345                  2179                  8084 
# S10_Hb_KDM_reanalysis S11_Hb_KDM_reanalysis S12_Hb_KDM_reanalysis 
# 6319                  8839                  4876 
# S3_Hb_KDM_reanalysis  S7_Hb_KDM_reanalysis  S8_Hb_KDM_reanalysis 
# 5050                  7802                  5317 
# S9_Hb_KDM_reanalysis 
# 7139 

message("Loaded seurat integrated.")


################## (2) load vector with valid barcodes by sample


# ## Scans arguments invoked from slurm job shell sh
# sample_tmp <- commandArgs(trailingOnly = TRUE)
# # For testing:
# # sample_tmp <- "4S_Hb_KDM_reanalysis, 4S_Hb_KDM"
# sample_data = unlist(strsplit(sample_tmp,","))
# Seurat_base_name <- trimws(sample_data[[2]])


## Read barcodes after remove outliers (PASS) and format cell-names 

unique(SeuratOBJ$orig.ident)[1]
strsplit(unique(SeuratOBJ$orig.ident)[1], split = "_")

csvDir_barcodes <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", "csv_files")
lst_bc <- list.files(csvDir_barcodes, pattern = "*_bc_PASS_isOutliers.csv")


all_bc_to_remove <- c()

for (bc_file in lst_bc) {
  
  message("Barcode list: ", bc_file)
  
  ## load barcodes that PASS outliers
  
  df_valid_barcodes_filtered <- read.csv(here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", "csv_files", bc_file))
  len_valid_bc <- length(df_valid_barcodes_filtered$x)
  v_valid_barcodes_filtered <- df_valid_barcodes_filtered$x
  head(v_valid_barcodes_filtered)
  message("Vector with ", len_valid_bc ," valid barcodes for sample `", bc_file,"` loaded")
  
  ## format barcodes to match seurat integrated barcode names
  
  ## short prefix
  prefixCell <- unlist(strsplit(bc_file, split = "_"))[1]
  barcodes_to_remove <- paste(prefixCell, "_", v_valid_barcodes_filtered, sep="")
  
  message("Cells to remove: ", length(barcodes_to_remove), " in sample ", prefixCell)
  
  all_bc_to_remove <- append(all_bc_to_remove, barcodes_to_remove) 
  
  ## long prefix
  #prefix_cell_name <- unlist(strsplit(bc_file, split = "_bc_PASS_isOutliers.csv"))[1]
  #prefix_cell_name <- paste0(prefix_cell_name,"_reanalysis")
  
}


################## (3) kept only cells with Outliers 

## remove cells from integrated seurat

message("Barcodes to remove: ", length(all_bc_to_remove))

SeuratObj_subset <- subset(SeuratOBJ, cells = all_bc_to_remove, invert = TRUE)

message("Resulting subset excluding cross-barcodes from Cell Ranger-count Seurat object: `", Seurat_base_name, "`")

print(table(SeuratObj_subset$orig.ident))
# 4S_Hb_KDM_reanalysis  5S_Hb_KDM_reanalysis  6S_Hb_KDM_reanalysis 
# 663                   256                   811 
# S10_Hb_KDM_reanalysis S11_Hb_KDM_reanalysis S12_Hb_KDM_reanalysis 
# 450                   625                   609 
# S3_Hb_KDM_reanalysis  S7_Hb_KDM_reanalysis  S8_Hb_KDM_reanalysis 
# 527                   393                   250 
# S9_Hb_KDM_reanalysis 
# 854 

#table(Idents(SeuratObj_subset))
rm("SeuratOBJ")


################## (4) Save Seurat object with ONLY Outliers to identify cell types with `01_Hb_celltypes_from_seurat_reanalyze.R` script


message("Final number of cells: ", sum(table(SeuratObj_subset$orig.ident)))

Seurat_base_name <- paste0(Seurat_base_name, "_subset_Outliers")

## Save new Seurat-subset 
rds_name <- paste0(Seurat_base_name, ".rds")
rds_name <- here(outputDir, rds_name)
saveRDS(SeuratObj_subset, file = rds_name) #seurat.norm_counts_Harmony_All_cellRanger_count_subset.rds

message("Saved Seurat subset data with ONLY Outliers!")


library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


