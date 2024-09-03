## Build csv file(s) to input in cellranger-ARC libraries (CSV files)
## It can handle samples with multiple gex or atac libraries, assuming they have the same base-name

# CellRanger-ARC library design
# fastqs,sample,library_type
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/raw-data/FASTQ/GEX/1HPC_S1_L003, 1G_HPC_KDM, Gene Expression
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/raw-data/FASTQ/ATAC/1HPC_S1_L003, 1A_HPC_KDM, Chromatin Accessibility
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/HPC_multiome_pilot/raw-data/FASTQ/ATAC/1HPC_S2_L003, 1A_HPC_KDM, Chromatin Accessibility
# ...
# Note can be more than one directory by side (atac/gex). In this case you need to adapt code


# Tracking multiome data from directories
library(stringr)
library(here)

here::here()

# Main FASTQ container
data_path <- here("raw-data", "FASTQ_2024_data_package2")
cellranger_library_path <- here("code", "cellranger")
# Read soft links for GEX Fastqs
snRNAseq_path <- here(data_path, "GEX")
dir_list <- list.files(path=snRNAseq_path, full.names=FALSE, recursive = FALSE)
dir_list

for (new_subdir in dir_list) {
  # new_subdir <-  "3C_Hb_KDM"
  # retrieve only unique IDs for each sample (includes R1, R2 and I1 and I2)
  dir_list <- list.files(path=here(snRNAseq_path, new_subdir), pattern="*R1", full.names=FALSE, recursive = FALSE)
  row_gexALL <- ""
  for (f in dir_list) {
    # f <- "3C_Hb_KDM_S2_L005_R1_001.fastq.gz"
    # Extract base-name following the cellranger-arc naming convention
    sub_dirID <- substring(f,1, regexpr("_L", f) + 4)
    # Append in rows
    fastqs_gex = here(snRNAseq_path, new_subdir, sub_dirID)
    row_gex <- paste0(fastqs_gex, ", ", new_subdir, ", Gene Expression\n" )
    row_gexALL <- paste(row_gexALL, row_gex)
  }
  cat(row_gexALL)

  ####### Append ATAC libraries
  snATACseq_path <- here(data_path, "ATAC")
  # Check corresponding ATAC directory exists. Assume same base-name is used for complementary atac side.
  new_subdir <- gsub("C", "A", new_subdir)
  fastqs_atac <- list.files(path=here(data_path, "ATAC", new_subdir), full.names=FALSE, recursive = FALSE)
  dir_list <- list.files(path=here(snATACseq_path, new_subdir), pattern="*R1", full.names=FALSE, recursive = FALSE)
  row_atacALL <- ""
  for (f in dir_list) {
    # f <- "3A_Hb_KDM_S3_L008_R1_001.fastq.gz"
    sub_dirID <- substring(f,1, regexpr("_L", f) + 4)
    # Append in rows
    fastqs_atac = here(snATACseq_path, new_subdir, sub_dirID)
    row_atac <- paste0(fastqs_atac, ", ", new_subdir, ", Chromatin Accessibility\n")
    row_atacALL <- paste(row_atacALL, row_atac)
  }
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

## Output sample:
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024_data_package2/GEX/9C_Hb_KDM/9C_Hb_KDM_S4_L005, 9C_Hb_KDM, Gene Expression
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024_data_package2/GEX/9C_Hb_KDM/9C_Hb_KDM_S4_L006, 9C_Hb_KDM, Gene Expression
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024_data_package2/GEX/9C_Hb_KDM/9C_Hb_KDM_S4_L007, 9C_Hb_KDM, Gene Expression
# /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/raw-data/FASTQ_2024_data_package2/ATAC/9A_Hb_KDM/9A_Hb_KDM_S5_L007, 9A_Hb_KDM, Chromatin Accessibility
# CellRanger-ARC library for : `9_Hb_KDM` saved on /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellranger//dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/code/cellranger/multiome_library_S9_Hb_KDM.csv



# 
# ## Parse the file directories
# for (new_subdir in dir_list) {
#   # values for testing
#   # new_subdir <- "4C_Hb_KDM"
#   
#   # # create header
#   # header <- 'fastqs,sample,library_type'
#   # # Append GEX rows
#   # fastqs_gex = here(snRNAseq_path, new_subdir)
#   # row_gex = paste0(fastqs_gex, ", ", new_subdir, ", Gene Expression")
#   # row_gex = paste0(header, '\n', row_gex) 
#   # row_gex
#   
#   # Append ATAC rows 
#   # We validate whether the corresponding ATAC directory exists
#   new_subdir <- gsub("C", "A", new_subdir)
#   fastqs_atac <- list.files(path=here(data_path, "ATAC", new_subdir), full.names=FALSE, recursive = FALSE)
#   
#   if (length(fastqs_atac) > 0) {
#     row_atac = paste0(here(snATACseq_path,new_subdir), ", ", new_subdir, ", Chromatin Accessibility")
#     row_all = paste0(row_gex, "\n", row_atac)
#     cat(row_all)
#     
#     ## Create and save library csv file
#     new_subdir <- gsub("A", "S", new_subdir)
#     library_name <- paste0("multiome_library_S",  new_subdir, '.csv')
#     library_name <- here(cellranger_library_path, library_name)
#     write.table(row_all, library_name, row.names = FALSE, col.names = FALSE, quote=FALSE)
#     message("CellRanger-ARC library for : `", new_subdir, "` saved on ", cellranger_library_path, '/')
#     
#   } else {
#     
#     ## This means only GEX dir available, but not ATAC dir, thus library is skipped
#     message("No ATAC directory `", new_subdir,"` analogous to GEX directory available. Library skipped!")
#     
#   }  
# }


