## Build csv file(s) to input in the cellranger-ARC pipeline

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
data_path <- here("raw-data", "FASTQ_2024")
cellranger_library_path <- here("code", "cellranger")
# Read soft links for GEX Fastqs
snRNAseq_path <- here(data_path, "GEX")
dir_list <- list.files(path=snRNAseq_path, full.names=FALSE, recursive = FALSE)
dir_list
# [1] "4C_Hb_KDM" "5C_Hb_KDM" "6C_Hb_KDM" "7C_Hb_KDM" "8C_Hb_KDM"
snATACseq_path <- here(data_path, "ATAC")
# dir_list_atac <- list.files(path=snATACseq_path, full.names=FALSE, recursive = FALSE)
# dir_list_atac
# [1] "3A_Hb_KDM" "4A_Hb_KDM" "5A_Hb_KDM" "6A_Hb_KDM"

## Parse the file directories
for (new_subdir in dir_list) {
  # values for testing
  # new_subdir <- "4C_Hb_KDM"
  
  # create header
  header <- 'fastqs,sample,library_type'
  # Append GEX rows
  fastqs_gex = here(snRNAseq_path, new_subdir)
  row_gex = paste0(fastqs_gex, ", ", new_subdir, ", Gene Expression")
  row_gex = paste0(header, '\n', row_gex) 
  row_gex
  
  # Append ATAC rows 
  # We validate whether the corresponding ATAC directory exists
  new_subdir <- gsub("C", "A", new_subdir)
  fastqs_atac <- list.files(path=here(data_path, "ATAC", new_subdir), full.names=FALSE, recursive = FALSE)
  
  if (length(fastqs_atac) > 0) {
    row_atac = paste0(here(snATACseq_path,new_subdir), ", ", new_subdir, ", Chromatin Accessibility")
    row_all = paste0(row_gex, "\n", row_atac)
    cat(row_all)
    
    ## Create and save library csv file
    new_subdir <- gsub("A", "S", new_subdir)
    library_name <- paste0("multiome_library_",  new_subdir, '.csv')
    library_name <- here(cellranger_library_path, library_name)
    #write.csv(row_all, library_name, row.names = FALSE, quote=FALSE)
    write.table(row_all, library_name, row.names = FALSE, col.names = FALSE, quote=FALSE)
    message("CellRanger-ARC library for : `", new_subdir, "` saved on ", cellranger_library_path, '/')
    
  } else {
    
    ## This means only GEX dir available, but not ATAC dir, thus library is skipped
    message("No ATAC directory `", new_subdir,"` analogous to GEX directory available. Library skipped!")
    
  }  
}


