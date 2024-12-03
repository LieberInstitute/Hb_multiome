########################################################################
## Prepare a Seurat RDS object WITHOUT OUTLIER cells
##  
## INPUT:
##      RDS (Seurat) with ONLY GEX or ATAC multiome cells with outliers
## OUPUT:
##      RDS (Seurat) with ONLY valid cells (good quality) detected (remove outliers)
## NOTE:
##      For +60k cells request 60G free-mem
##
## Authors. CSC 
## Date. Dec 2024
########################################################################

library("Seurat")
library("Signac") 
library("here")
library("tidyr")
library("stringr")

here::here()

# Check/create directories
inputDir_fulldataset <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze")
inputDir_outliers <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze_outliers") 
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
if (length(list.files(inputDir_outliers, pattern = Seurat_base_name)==1)) {
  message("Processing ", Seurat_base_name)
} else {
  message("Input seurat object missed!")
  stop()
}

message("Loading atypicals for multiome RNA `", Seurat_base_name, "`")

## Load Seurat Integrated with cluster information
SeuratOBJ <- readRDS(here(inputDir_outliers, Seurat_base_name))
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
if (length(list.files(inputDir_outliers, pattern = Seurat_base_name)==1)) {
  message("Processing ", Seurat_base_name)
} else {
  message("Input seurat object missed!")
  stop()
}

message("Loading atypicals for multiome ATAC `", Seurat_base_name, "`")

## Load Seurat Integrated with cluster information
SeuratOBJ <- readRDS(here(inputDir_outliers, Seurat_base_name))
## verification
table(SeuratOBJ$orig.ident)
# 4S_Hb_KDM_reanalysis  5S_Hb_KDM_reanalysis  6S_Hb_KDM_reanalysis 
# 539                    35                   671 
# S10_Hb_KDM_reanalysis S11_Hb_KDM_reanalysis S12_Hb_KDM_reanalysis 
# 1                    26                    98 
# S3_Hb_KDM_reanalysis  S8_Hb_KDM_reanalysis 
# 12                    84 

all_atac_bc_to_keep <- Cells(SeuratOBJ)
rm("SeuratOBJ")
message(length(all_atac_bc_to_keep), " cells detected as atypicals in ATAC")
# 1466 cells detected as atypicals in RNA.


## get barcodes from both RNA and ATAC outliers
all_rna_atac_bc_to_keep <- union(all_rna_bc_to_keep, all_atac_bc_to_keep)

message(length(all_rna_atac_bc_to_keep), " common barcodes to remove")
#7248 common barcodes to remove



################## (3) remove atypical from cellrangerARC-reanalyze dataset

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts'}
Seurat_base_name <- paste0(Seurat_base_name, '_Harmony_All')

## Validate seurat exists
if (length(list.files(inputDir_fulldataset, pattern = Seurat_base_name)==1)) {
  message("Processing ", Seurat_base_name)
} else {
  message("Input seurat object missed!")
  stop()
}

message("Removing atypical on ", Seurat_base_name)

## Load Seurat Integrated with cluster information
SeuratOBJ <- readRDS(here(inputDir_fulldataset, paste0(Seurat_base_name, ".rds")))
table(SeuratOBJ$orig.ident)
# 4S_Hb_KDM_reanalysis  5S_Hb_KDM_reanalysis  6S_Hb_KDM_reanalysis 
# 7345                  2179                  8084 
# S10_Hb_KDM_reanalysis S11_Hb_KDM_reanalysis S12_Hb_KDM_reanalysis 
# 6319                  8839                  4876 
# S3_Hb_KDM_reanalysis  S7_Hb_KDM_reanalysis  S8_Hb_KDM_reanalysis 
# 5050                  7802                  5317 
# S9_Hb_KDM_reanalysis 
# 7139 

message("Removing outliers from ", Seurat_base_name)

SeuratObj_subset <- subset(SeuratOBJ, cells = all_rna_atac_bc_to_keep, invert = TRUE)
rm("SeuratOBJ")

x <- as.data.frame(table(SeuratObj_subset$orig.ident))[2]
rownames(x) <- NULL
print(x)
# 4S_Hb_KDM_reanalysis  5S_Hb_KDM_reanalysis  6S_Hb_KDM_reanalysis 
# 6364                  1990                  7059 
# S10_Hb_KDM_reanalysis S11_Hb_KDM_reanalysis S12_Hb_KDM_reanalysis 
# 5914                  8054                  4045 
# S3_Hb_KDM_reanalysis  S7_Hb_KDM_reanalysis  S8_Hb_KDM_reanalysis 
# 4375                  7207                  4735 
# S9_Hb_KDM_reanalysis 
# 5959 

message(sum(x), " cells retained")
# 55702 cells retained


################## (4) Rename sample IDs and cell names before clustering (to easy sorting tasks)


## Rename samples IDs
unique(SeuratObj_subset@meta.data$orig.ident)
SeuratObj_subset$orig.ident <- sprintf("S%02d_Hb_r", readr::parse_number(SeuratObj_subset$orig.ident))
unique(SeuratObj_subset@meta.data$orig.ident)
# [1] "S04_Hb_r" "S05_Hb_r" "S06_Hb_r" "S10_Hb_r" "S11_Hb_r" "S12_Hb_r"
# [7] "S03_Hb_r" "S07_Hb_r" "S08_Hb_r" "S09_Hb_r"

## Rename cell names
Cells(SeuratObj_subset)[1:5]
# [1] "4S_AAACAGCCAGAATGAC-1" "4S_AAACAGCCAGCAAGGC-1" "4S_AAACATGCACCTGGTG-1"
# [4] "4S_AAACATGCAGGATGGC-1" "4S_AAACATGCAGTAATAG-1"
tail(Cells(SeuratObj_subset), n=5)
# [1] "S9_TTTGTGTTCCGTGACA-1" "S9_TTTGTGTTCCGTTATT-1" "S9_TTTGTGTTCGTTAGCG-1"
# [4] "S9_TTTGTTGGTCATGCAA-1" "S9_TTTGTTGGTTGTTCAC-1"
## Prepare vector with new cell names
v_new_cell_names <-paste0(sprintf("S%02d", readr::parse_number(Cells(SeuratObj_subset))), "_", 
                     str_extract(Cells(SeuratObj_subset), regex("[ACTG]*-1")))
length(v_new_cell_names)
## rename cells 
SeuratObj_subset <- RenameCells(SeuratObj_subset, new.names = v_new_cell_names)
Cells(SeuratObj_subset)[1:5]
tail(Cells(SeuratObj_subset), n=5)

message(length(Cells(SeuratObj_subset)), " cells retained and formatted for further analysis.")



################## (4) Find DEG and save Seurat with ONLY ATAC OUTLIER cells (barcodes) to identify cell types later

## Find DEG in the integrated Seurat for ALL clusters (BEFORE pseudobulk)
#table(SeuratOBJ[["seurat_clusters"]])
all.markers <- FindAllMarkers(object = SeuratObj_subset)
#head(all.markers, n=3)

cvs_file <- paste0(Seurat_base_name, "QCed_markers_RNA_ATAC.csv") 
cvs_file <- here(cvsDir, cvs_file)
write.csv(all.markers, cvs_file)

message(" FindAllMarkers done!")

## Save Seurat without outlier cells

# ## Format sample IDs and cell-name IDs
# 
# Cells(SeuratObj_subset)[1:10]
# extract_numeric(Cells(SeuratObj_subset)[1:10])
# df_mdT$orig.ident  <- sprintf("S%02d_Hb_r", extract_numneric(df_mdT$orig.ident))
# df_mdT$orig.ident  <- sprintf("S%02d_Hb_r", extract_numneric(df_mdT$orig.ident))

Seurat_base_name <- paste0(Seurat_base_name, "_ARC_reanalize_QCed")

## Save new Seurat-subset
rds_name <- paste0(Seurat_base_name, ".rds")
rds_name <- here(outputDir, rds_name)
saveRDS(SeuratObj_subset, file = rds_name)

message("Saved Seurat subset data QCed!")





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


