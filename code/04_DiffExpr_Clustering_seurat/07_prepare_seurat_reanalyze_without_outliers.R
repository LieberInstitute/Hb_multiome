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

# Check/create directories
not_filtered_inputDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze")
inputDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze_outliers") 
outputDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze")
cvsDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze_outliers", "cvs_files_markers")

## Check directories
if (!dir.exists(outputDir)) {dir.create(outputDir)}
if (!dir.exists(cvsDir)) {dir.create(cvsDir)}


#############################           Initials        ################################

## Set count-mtx type and integration model (CCA or Harmony)

count_mtx_type <- 'norm_counts' 
Seurat_reduction <- 'Harmony' 
minCells <- 1

################## (1) load vector with valid barcodes by sample for RNA assays

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts'}
Seurat_base_name <- paste0(Seurat_base_name, "_Harmony_All_GEX_subset_Outliers.rds")

## Validate seurat exists
if (length(list.files(inputDir, pattern = Seurat_base_name)==1)) {
  message("Processing ", Seurat_base_name)
} else {
  message("Input seurat object missed!")
  stop()
}

message("Loading atypicals for multiome RNA `", Seurat_base_name, "`")

## Load Seurat Integrated with cluster information
SeuratOBJ <- readRDS(here(inputDir, Seurat_base_name))
## verification
table(SeuratOBJ$orig.ident)
# 4S_Hb_KDM_reanalysis  5S_Hb_KDM_reanalysis  6S_Hb_KDM_reanalysis 
# 564                   158                   449 
# S10_Hb_KDM_reanalysis S11_Hb_KDM_reanalysis S12_Hb_KDM_reanalysis 
# 405                   762                   751 
# S3_Hb_KDM_reanalysis  S7_Hb_KDM_reanalysis  S8_Hb_KDM_reanalysis 
# 664                   595                   515 
# S9_Hb_KDM_reanalysis 
# 1180 

all_rna_bc_to_keep <- Cells(SeuratOBJ)
rm("SeuratOBJ")
message(length(all_rna_bc_to_keep), " cells detected as atypicals in RNA.")
# 6043 cells detected as atypicals in RNA.



################## (2) load vector with valid barcodes for ATAC assays

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts'}
Seurat_base_name <- paste0(Seurat_base_name, "_Harmony_All_ATAC_subset_Outliers.rds")

## Validate seurat exists
if (length(list.files(inputDir, pattern = Seurat_base_name)==1)) {
  message("Processing ", Seurat_base_name)
} else {
  message("Input seurat object missed!")
  stop()
}

message("Loading atypicals for multiome ATAC `", Seurat_base_name, "`")

## Load Seurat Integrated with cluster information
SeuratOBJ <- readRDS(here(inputDir, Seurat_base_name))
## verification
# table(SeuratOBJ$orig.ident)
# 4S_Hb_KDM_reanalysis  5S_Hb_KDM_reanalysis  6S_Hb_KDM_reanalysis 
# 539                    35                   671 
# S10_Hb_KDM_reanalysis S11_Hb_KDM_reanalysis S12_Hb_KDM_reanalysis 
# 1                    26                    98 
# S3_Hb_KDM_reanalysis  S8_Hb_KDM_reanalysis 
# 12                    84 

all_atac_bc_to_keep <- Cells(SeuratOBJ)
rm("SeuratOBJ")
message(length(all_atac_bc_to_keep), " cells detected as atypicals in RNA.")
# 1466 cells detected as atypicals in RNA.

length(all_atac_bc_to_keep)
# [1] 61484 ATAC


all_rna_atac_bc_to_keep <- union(all_rna_bc_to_keep, all_atac_bc_to_keep)

message(length(all_rna_atac_bc_to_keep), " common barcodes to remove")
#7248 common barcodes to remove



################## (3) remove atypical from dataset

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts'}
Seurat_base_name <- paste0(Seurat_base_name, '_Harmony_All')

## Validate seurat exists
if (length(list.files(not_filtered_inputDir, pattern = Seurat_base_name)==1)) {
  message("Processing ", Seurat_base_name)
} else {
  message("Input seurat object missed!")
  stop()
}

message("Removing atypical on ", Seurat_base_name)

## Load Seurat Integrated with cluster information
SeuratOBJ <- readRDS(here(not_filtered_inputDir, paste0(Seurat_base_name, ".rds")))
table(SeuratOBJ$orig.ident)
# 4S_Hb_KDM_reanalysis  5S_Hb_KDM_reanalysis  6S_Hb_KDM_reanalysis 
# 7345                  2179                  8084 
# S10_Hb_KDM_reanalysis S11_Hb_KDM_reanalysis S12_Hb_KDM_reanalysis 
# 6319                  8839                  4876 
# S3_Hb_KDM_reanalysis  S7_Hb_KDM_reanalysis  S8_Hb_KDM_reanalysis 
# 5050                  7802                  5317 
# S9_Hb_KDM_reanalysis 
# 7139 

SeuratObj_subset <- subset(SeuratOBJ, cells = all_rna_atac_bc_to_keep, invert = TRUE)
rm("SeuratOBJ")

table(SeuratObj_subset$orig.ident)

## Save Seurat with outlier cells

Seurat_base_name <- paste0(Seurat_base_name, "_ARC_reanalize_qced")

## Save new Seurat-subset
rds_name <- paste0(Seurat_base_name, ".rds")
rds_name <- here(outputDir, rds_name)
saveRDS(SeuratObj_subset, file = rds_name)

message("Saved Seurat subset data QCed!")




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


