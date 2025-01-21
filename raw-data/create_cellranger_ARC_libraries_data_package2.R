## Build csv file(s) to input in cellranger-ARC libraries (CSV files)
## It can handle samples with multiple gex or atac libraries, assuming they have the same base-name

# CellRanger-ARC library naming-convention requirements for cellranger-arc v2+ (one line by library)
# fastqs,sample,library_type
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/raw-data/FASTQ/GEX/1HPC_S1_L003, 1G_HPC_KDM, Gene Expression
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/raw-data/FASTQ/ATAC/1HPC_S1_L003, 1A_HPC_KDM, Chromatin Accessibility
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/raw-data/FASTQ/ATAC/1HPC_S2_L003, 1A_HPC_KDM, Chromatin Accessibility
# ...
# Note can be more than one directory by side (atac/gex). In this case you need to adapt code


# Tracking multiple fastq libraries from raw-data directories
library(stringr)
library(here)

here::here()

# Main FASTQ container
data_path <- here("raw-data", "FASTQ_2024_data_package2")
cellranger_library_path <- here("code", "cellranger")
# Read soft links for GEX Fastqs
snRNAseq_path <- here(data_path, "GEX")
dir_list <- list.files(path=snRNAseq_path, full.names=FALSE, recursive = FALSE)

for (new_subdir in dir_list) {
  
  # new_subdir <-  "9C_Hb_KDM"
  row_gexALL <- ""
  
  ### # for create one line for gex library
  # Append all in one row
  fastqs_gex = here(snRNAseq_path, new_subdir)
  row_gexALL <- paste0(fastqs_gex, ", ", new_subdir, ", Gene Expression\n" )
  cat(row_gexALL)

  ####### Append ATAC libraries
  snATACseq_path <- here(data_path, "ATAC")
  # Check corresponding ATAC directory exists. Assume same base-name is used for complementary atac side.
  new_subdir <- gsub("C", "A", new_subdir)
  row_atacALL <- ""
  
  # Append all in one row
  fastqs_atac = here(snATACseq_path, new_subdir)
  row_atac <- paste0(fastqs_atac, ", ", new_subdir, ", Chromatin Accessibility\n")
  row_atacALL <- paste0(row_atacALL, row_atac)
  cat(row_atacALL)
  
  # create body csv text
  header <- 'fastqs,sample,library_type\n'
  body_cvs <- paste0(header, row_gexALL, row_atacALL) 
  cat(body_cvs)  
  
  ## Create and save library csv file
  new_subdir <- gsub("A", "", new_subdir)
  library_name <- paste0("multiome_library_S",  new_subdir, '.csv')
  library_name <- here(cellranger_library_path, library_name)
  write.table(body_cvs, library_name, row.names = FALSE, col.names = FALSE, quote=FALSE)
  message("CellRanger-ARC library for : `", new_subdir, "` saved on ", cellranger_library_path, '/', library_name)
  
  
}

## Output 1. example for one sample with multiple gex libraries:
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024_data_package2/GEX/9C_Hb_KDM/9C_Hb_KDM_S4_L005, 9C_Hb_KDM, Gene Expression
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024_data_package2/GEX/9C_Hb_KDM/9C_Hb_KDM_S4_L006, 9C_Hb_KDM, Gene Expression
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024_data_package2/GEX/9C_Hb_KDM/9C_Hb_KDM_S4_L007, 9C_Hb_KDM, Gene Expression
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024_data_package2/ATAC/9A_Hb_KDM/9A_Hb_KDM_S5_L007, 9A_Hb_KDM, Chromatin Accessibility
# CellRanger-ARC library for : `9_Hb_KDM` saved on /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellranger//dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellranger/multiome_library_S9_Hb_KDM.csv

## Output 2. example for one sample with all libraries in the same directory:
# fastqs,sample,library_type
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024_data_package2/GEX/3C_Hb_KDM, 3C_Hb_KDM, Gene Expression
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024_data_package2/ATAC/3A_Hb_KDM, 3A_Hb_KDM, Chromatin Accessibility


## This chunk is ONLY for sample S7 and S8 (they arrived apart)

dir_list <- c("7C_Hb_KDM", "8C_Hb_KDM")
for (new_subdir in dir_list) {
  row_gexALL <- ""
  # Append all in one row
  fastqs_gex = here(snRNAseq_path, new_subdir)
  row_gexALL <- paste0(fastqs_gex, ", ", new_subdir, ", Gene Expression\n" )
  cat(row_gexALL)
  
  ####### Append ATAC libraries
  snATACseq_path <- here(data_path, "ATAC")
  # Check corresponding ATAC directory exists. Assume same base-name is used for complementary atac side.
  new_subdir <- gsub("C", "A", new_subdir)
  row_atacALL <- ""
  
  # Append all in one row
  fastqs_atac = here(snATACseq_path, new_subdir)
  row_atac <- paste0(fastqs_atac, ", ", new_subdir, ", Chromatin Accessibility\n")
  row_atacALL <- paste0(row_atacALL, row_atac)
  cat(row_atacALL)
  
  # create body csv text
  header <- 'fastqs,sample,library_type\n'
  body_cvs <- paste0(header, row_gexALL, row_atacALL) 
  cat(body_cvs)  
  
  ## Create and save library csv file
  new_subdir <- gsub("A", "", new_subdir)
  library_name <- paste0("multiome_library_S",  new_subdir, '.csv')
  library_name <- here(cellranger_library_path, library_name)
  write.table(body_cvs, library_name, row.names = FALSE, col.names = FALSE, quote=FALSE)
  message("CellRanger-ARC library for : `", new_subdir, "` saved on ", cellranger_library_path, '/', library_name)
  
}


