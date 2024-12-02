########################################################################
## Prepare a Seurat object WITHOUT OUTLIERS to identify cell-types for both GEX and ATAC muliome dataset
##  
## INPUT:
##      CSV files with barcodes that passed outliers calculated on both GEX and ATAC data
## 
## OUPUT:
##      Seurat with valid cells (good quality) detected
## NOTE:
##      For +60k cells request 60G free-mem
##
## Authors. CSC 
## Date. Dec 2024
########################################################################

library("Seurat")
library("Signac") 
library("here")

here::here()


################## (1) Load integrated Seurat with Harmony correction


# Check/create directories
inputDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze")
outputDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze_outliers") 
cvsDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze_outliers", "cvs_files_markers")

## Check directories
if (!dir.exists(outputDir)) {dir.create(outputDir)}
if (!dir.exists(cvsDir)) {dir.create(cvsDir)}


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
SeuratOBJ <- readRDS(here(inputDir, paste0(Seurat_base_name, ".rds")))
## verification
table(SeuratOBJ$orig.ident)
# 4S_Hb_KDM_reanalysis  5S_Hb_KDM_reanalysis  6S_Hb_KDM_reanalysis 
# 7345                  2179                  8084 
# S10_Hb_KDM_reanalysis S11_Hb_KDM_reanalysis S12_Hb_KDM_reanalysis 
# 6319                  8839                  4876 
# S3_Hb_KDM_reanalysis  S7_Hb_KDM_reanalysis  S8_Hb_KDM_reanalysis 
# 5050                  7802                  5317 
# S9_Hb_KDM_reanalysis 
# 7139 

message("Loaded seurat integrated with ", length(Cells(SeuratOBJ)), " cells.")
# Loaded seurat integrated with 62950 cells




################## (2) load vector with valid barcodes by sample ONLY for GEX assays

## Read barcodes and build a Seurat to keep only barcodes detected as outliers. Format cell-names too. 

#unique(SeuratOBJ$orig.ident)[1]

csvDir_barcodes <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", "csv_files")
lst_bc <- list.files(csvDir_barcodes, pattern = "*_bc_PASS_GEX_isOutliers.csv")
# [1] "4S_Hb_KDM_bc_PASS_GEX_isOutliers.csv" 
# [2] "5S_Hb_KDM_bc_PASS_GEX_isOutliers.csv" 
# [3] "6S_Hb_KDM_bc_PASS_GEX_isOutliers.csv" 
# [4] "S10_Hb_KDM_bc_PASS_GEX_isOutliers.csv"
# [5] "S11_Hb_KDM_bc_PASS_GEX_isOutliers.csv"
# [6] "S12_Hb_KDM_bc_PASS_GEX_isOutliers.csv"
# [7] "S3_Hb_KDM_bc_PASS_GEX_isOutliers.csv" 
# [8] "S7_Hb_KDM_bc_PASS_GEX_isOutliers.csv" 
# [9] "S8_Hb_KDM_bc_PASS_GEX_isOutliers.csv" 
# [10] "S9_Hb_KDM_bc_PASS_GEX_isOutliers.csv" 

all_rna_bc_to_keep <- c()

## Parse csv files containing valid barcodes

for (bc_file in lst_bc) {
  # testing: bc_file <- lst_bc[1]
  message("Valid barcodes list: ", bc_file)
  
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
  
  message("Cells to keep on multiome GEX side: ", length(barcodes_to_remove), " in sample ", prefixCell)
  
  all_rna_bc_to_keep <- append(all_rna_bc_to_keep, barcodes_to_remove) 
  
}

length(all_rna_bc_to_keep)


################## (3) kept only cells with Outliers 

## remove cells from integrated seurat

message("Valid barcodes on RNA side: ", length(all_rna_bc_to_keep))
# Valid barcodes on RNA side: 56907

# SeuratObj_subset <- subset(SeuratOBJ, cells = all_bc_to_keep)

# message("CellRangerARC-reanalyze GEX Seurat with ONLY outliers done!")

# print(table(SeuratObj_subset$orig.ident))

# message(length(Cells(SeuratObj_subset)), " valid barcodes for RNA side")

#table(Idents(SeuratObj_subset))
#rm("SeuratOBJ")



# ################## (4) Find DEG and save Seurat with ONLY OUTLIER cells (barcodes) to identify cell types later
# 
# ## Find DEG in the integrated Seurat for ALL clusters (BEFORE pseudobulk)
# #table(SeuratOBJ[["seurat_clusters"]])
# all.markers <- FindAllMarkers(object = SeuratObj_subset)
# #head(all.markers, n=3)
# 
# # cvs_file <- paste0(Seurat_base_name, '_', integration_model, '_Allmarkers.csv')
# cvs_file <- paste0(Seurat_base_name, "markers_GEX.csv") 
# #seurat.norm_counts_Harmony_Allmarkers_GEX.csv
# cvs_file <- here(cvsDir, cvs_file)
# write.csv(all.markers, cvs_file)
# 
# message(" FindAllMarkers done!")
# 
# ## Save Seurat with outlier cells
# 
# Seurat_base_name <- paste0(Seurat_base_name, "_GEX_subset_Outliers")
# 
# ## Save new Seurat-subset 
# rds_name <- paste0(Seurat_base_name, ".rds")
# rds_name <- here(outputDir, rds_name)
# saveRDS(SeuratObj_subset, file = rds_name) 
# 
# message("Saved Seurat subset data with ONLY Outliers!")



################## (5) load vector with valid barcodes by sample ONLY for ATAC assays

## Read barcodes and build a Seurat to keep only barcodes detected as outliers. Format cell-names too. 

#unique(SeuratOBJ$orig.ident)[1]

csvDir_barcodes <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", "csv_files")
lst_bc <- list.files(csvDir_barcodes, pattern = "*_bc_PASS_ATAC_isOutliers.csv")
# [1] "4S_Hb_KDM_bc_PASS_ATAC_isOutliers.csv" 
# [2] "5S_Hb_KDM_bc_PASS_ATAC_isOutliers.csv" 
# [3] "6S_Hb_KDM_bc_PASS_ATAC_isOutliers.csv" 
# [4] "S10_Hb_KDM_bc_PASS_ATAC_isOutliers.csv"
# [5] "S11_Hb_KDM_bc_PASS_ATAC_isOutliers.csv"
# [6] "S12_Hb_KDM_bc_PASS_ATAC_isOutliers.csv"
# [7] "S3_Hb_KDM_bc_PASS_ATAC_isOutliers.csv" 
# [8] "S7_Hb_KDM_bc_PASS_ATAC_isOutliers.csv" 
# [9] "S8_Hb_KDM_bc_PASS_ATAC_isOutliers.csv" 
# [10] "S9_Hb_KDM_bc_PASS_ATAC_isOutliers.csv" 

all_bc_to_keep <- c()

## Track total number of valid cells
names(table(SeuratOBJ$orig.ident))
as.list(as.vector(table(SeuratOBJ$orig.ident)))

## build the vector of cells to remove to build the ATAC isOutliers dataset 

for (bc_file in lst_bc) {
  # bc_file <- lst_bc[1]
  message("Barcode list: ", bc_file)
  
  ## load barcodes that PASS outliers ATAC
  
  df_valid_barcodes_filtered <- read.csv(here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", "csv_files", bc_file))
  len_valid_bc <- length(df_valid_barcodes_filtered$x)
  v_valid_barcodes_filtered <- df_valid_barcodes_filtered$x
  head(v_valid_barcodes_filtered)
  message("Vector with ", len_valid_bc ," valid barcodes for sample `", bc_file,"` loaded")
  
  ## format barcodes to match seurat integrated barcode names
  
  ## short prefix
  prefixCell <- unlist(strsplit(bc_file, split = "_"))[1]
  barcodes_to_remove <- paste(prefixCell, "_", v_valid_barcodes_filtered, sep="")
  
  message("Cells to remove on multiome GEX side: ", length(barcodes_to_remove), " in sample ", prefixCell)
  
  all_bc_to_keep <- append(all_bc_to_keep, barcodes_to_remove)
  
}

length(all_bc_to_keep)
# [1] 61484 ATAC

## remove cells from integrated seurat

message("Barcodes to remove on ATAC side: ", length(all_bc_to_keep))

table(SeuratOBJ$orig.ident)
SeuratObj_subset_atac <- subset(SeuratOBJ, cells = all_bc_to_keep, invert = TRUE)
table(SeuratObj_subset_atac$orig.ident)

message("CellRangerARC-reanalyze ATAC Seurat with ONLY outliers done!")



message(length(Cells(SeuratObj_subset)), " outliers found on ATAC")

#table(Idents(SeuratObj_subset))
#rm("SeuratOBJ")


################## (4) Find DEG and save Seurat with ONLY ATAC OUTLIER cells (barcodes) to identify cell types later

## Find DEG in the integrated Seurat for ALL clusters (BEFORE pseudobulk)
#table(SeuratOBJ[["seurat_clusters"]])
all.markers <- FindAllMarkers(object = SeuratObj_subset)
#head(all.markers, n=3)

cvs_file <- paste0(Seurat_base_name, "markers_ATAC.csv") 
#seurat.norm_counts_Harmony_Allmarkers_GEX.csv
cvs_file <- here(cvsDir, cvs_file)
write.csv(all.markers, cvs_file)

message(" FindAllMarkers done!")

## Save Seurat with outlier cells

Seurat_base_name <- paste0(Seurat_base_name, "_ATAC_subset_Outliers")

## Save new Seurat-subset 
rds_name <- paste0(Seurat_base_name, ".rds")
rds_name <- here(outputDir, rds_name)
saveRDS(SeuratObj_subset, file = rds_name) 

message("Saved Seurat subset data with ONLY Outliers!")





# # slurm script reproducibility
# library("slurmjobs")
# job_single(
#   name = "06_prepare_seurat_reanalyze_with_outliers", memory = "30G", cores = 2, create_shell = TRUE
# )

library("sessioninfo")
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()


