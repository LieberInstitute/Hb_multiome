#########################################################################
##
## Makeup DEG CVS files in a nice format
##
## Input:  All the files in the corresponding folder containing the CSV files
## Output: Files with new format suitable to share or add as supplemental material
## 
## Authors. CSC 
## Date: Feb 2025
########################################################################

library("sessioninfo")
library("tidyverse")
library("here")

here::here()

## Preparing directories

# path before cell annotations
# inputCVS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method", "cvs_files")
# outputCVS_Dir <- here("processed-data", "05_Clustering_ARCr", "01_clustering_std_method", "cvs_files_formatted")

inputCVS_Dir <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze_v3", "cvs_files_markers")
outputCVS_Dir <- here("processed-data", "05_Clustering_ARCr", "02_Hb_celltypes_from_seurat_reanalyze_v3", "cvs_files_formatted")

if (!dir.exists(outputCVS_Dir)) { dir.create(outputCVS_Dir) }

## Previous WNN inspected
lst_files <- list.files(inputCVS_Dir, pattern = "*_top50.csv$", full.names = TRUE, include.dirs = FALSE) 
basename(lst_files)
# [1] "WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r1_cellTypes_integrated_top50.csv"
# [2] "WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_cellTypes_integrated_top50.csv"
# [3] "WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r1_cellTypes_integrated_top50.csv"
# [4] "WNN_rnaHarm_atacHarm_k40_C.leiden_lsi_r2_cellTypes_integrated_top50.csv"

# Now, I only redo the annotations nice for our target WNN result
lst_files <- lst_files[basename(lst_files) == "WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_cellTypes_integrated_top50.csv"]
basename(lst_files)


# We do not need this any more
# lst_idx_to_remove = list()
# length(lst_files)
# for (f in seq_along(lst_files)) {
#   # f = 2
#   if (file.size(lst_files[f]) <= 3) { 
#     print(basename(lst_files[f])) 
#     lst_idx_to_remove <- append(lst_idx_to_remove, f)
#   }
# }
# ## remove empty files
# if (!length(lst_idx_to_remove) == 0) {
#   unlist(lst_idx_to_remove)
#   lst_files <- lst_files[-unlist(lst_idx_to_remove)]
#   length(lst_files)
# }


## rename and re-arrange columns for supplemental file (paper)

for (f in lst_files) {
  # f= lst_files[1]  
  message("renaming ", basename(f))
  df <- read.csv(f, header = T, sep = ",")
  #colnames(df)
  names(df)[names(df) == "p_val"] <- "p-value"
  names(df)[names(df) == "avg_log2FC"]  = "Log-fold-change (vs all clusters)"
  names(df)[names(df) == "pct.1"] <- "Fraction of cell types expressing"
  names(df)[names(df) == "pct.2"]  = "Fraction of all other cells expressing"
  names(df)[names(df) == "p_val_adj"] <- "FDR adjusted p-value"
  names(df)[names(df) == "cluster"]  = "Cell type ID"
  names(df)[names(df) == "gene"] <- "Gene"
  names(df)[names(df) == "cell_type"] <- "Marker-type"
  #colnames(df)

  # rearrange columns and save new formatted file
  df <- df[, c("Cell type ID", "Gene", "p-value", "Log-fold-change (vs all clusters)", 
               "Fraction of cell types expressing", "Fraction of all other cells expressing", "FDR adjusted p-value", "Marker-type")]
  write.csv(df, here(outputCVS_Dir, basename(f)), row.names = F)
  
}


##  slurm script reproducibility
# library("slurmjobs")
# job_single(
#   name = "10_format_DEG_csv_files", memory = "20G", cores = 1, create_shell = TRUE, partition = "katun", logdir = "logs", create_logdir = FALSE
# )

## Reproducibility information

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
