library(here)

################## (1) Load integrated Seurat with Harmony correction


################## (2) load vector with valid barcodes by sample

here::here()

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



