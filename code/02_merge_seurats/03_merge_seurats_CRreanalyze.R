########################################################################
## Merge Seurat objects contained in an specific directory
##
## Authors. CSC
## Date. Sep 25, 2024
## Last.md: XXX
##
## Input: Seurat RDS Object generated with 05_merge_seurats_job_loop_reanalyze.R
## Output:  New Seurat objects merged 
##
## NOTES: 80G free mem recommended for 60 to 80 thousand cells
## For slurm env: runsrun --x11 --pty --partition=interactive bash
########################################################################

library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(here)

here::here()

## processed_data and plots directories 
processedDir_out <- here("processed-data", "02_merge_seurats", "cellrangerARC_reanalyze")
if (!dir.exists(processedDir_out)) { dir.create(processedDir_out) }


########################    Initials ########################  

## Counts for GEX assay available in Seurat: raw and normalized data
##        count_mtx_type_label <- 'data_counts'      
##        count_mtx_type_label <- 'norm_counts' 

args = commandArgs(trailingOnly=TRUE)
## read count mtx type:
##      all_type_mtx=(data_counts norm_counts)

count_mtx_type_label <- args[2]
## for testing:  count_mtx_type_label <- "norm_counts" 

## compose rds base name for merge data
if (count_mtx_type_label=='data_counts') { 
  s_sample <- 'seurat.combined.data_counts' 
  rna_layer <- 'counts'
} else if (count_mtx_type_label=='norm_counts') { 
  s_sample <- 'seurat.combined.norm_counts' 
  rna_layer <- 'data'
} else {
  message("Assay type not provided!")
  stop()
}

message('Merging  `', count_mtx_type_label, '` seurat assays.')
message('\nSuffix name for merged Seurat is `', s_sample, "`")

## read directory with Seurat objects
processedDir_in <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze")
all_rds <- here(processedDir_in, list.files(processedDir_in, pattern="_reanalysis.rds"))

if (length(all_rds) < 1) { stop('\nOnly one Seurat available.') }

message("Seurats to combine: ")
print(basename(all_rds))


### Prepare list of Seurat(s) to merge

seurat_lst <- list()
seurat_name_lst <- list()

for (rds_path in all_rds) {
  
  print(rds_path)
  SeuratOBJ <-readRDS(rds_path)
  print(SeuratOBJ)
  DefaultAssay(SeuratOBJ) <- "RNA"
  
  if (count_mtx_type_label=='norm_counts') { SeuratOBJ <- NormalizeData(SeuratOBJ) }
  
  if ( length(seurat_lst)>0 ) { seurat_lst <- append(seurat_lst, SeuratOBJ) } else { seurat_lst <- SeuratOBJ }
  # Get first word for the Seurat name
  s <- strsplit(levels(SeuratOBJ$orig.ident)[1], split = "_")[[1]][1]
  if ( length(seurat_name_lst)>0 ) { seurat_name_lst <- append(seurat_name_lst, s) } else { seurat_name_lst <- s }
  
}

# print(seurat_lst)
# # testing
# Cells(SeuratOBJ)[1:10]
# Features(SeuratOBJ)
# nrow(SeuratOBJ)

# some validations, 'data' slot should exists if NormalizedData was ran
# Layers(SeuratOBJ)
# max(SeuratOBJ[["RNA"]]$counts) #1255
# max(SeuratOBJ[["RNA"]]$data) #8.2943


message('\nYou combined ', length(seurat_lst), ' Seurat objects: ', sapply(seurat_name_lst, function(i) paste0(i, ', ')) )



## Merge the Seurat objects contained in the list according with the `count type`
## NOTE: By default, merge() will combine the Seurat objects based on the raw count matrices, erasing any previously normalized and scaled data matrices. If you want to merge the normalized data matrices as well as the raw count matrices, simply pass merge.data = TRUE. This should be done if the same normalization approach was applied to all objects.

if (length(seurat_lst) > 1) {
  
  # assign first Seurat to merge the objects  
  SeuratOBJ <- seurat_lst[[1]]
  l <- length(seurat_lst)
  # assign the x Seurats to merge
  #if (l>2) { SeuratOBJx <- (seurat_lst)[[2:l]] } else { SeuratOBJx <- seurat_lst[[l]] }
  if (l>2) { SeuratOBJx <- (seurat_lst)[2:l] } else { SeuratOBJx <- seurat_lst[l] }
  
  # By default, merge() will combine the Seurat objects based on the raw count matrices
  SeuratOBJ.combined <- merge(SeuratOBJ, y = c(SeuratOBJx), 
                              add.cell.ids = c(seurat_name_lst), 
                              project = "Habenula",
                              merge.data = TRUE)     #  merge the normalized and raw count 
  ##pbmc.big <- merge(pbmc3k, y = c(pbmc4k, pbmc8k), add.cell.ids = c("3K", "4K", "8K"), project = "PBMC15K")
  
} else {
  message('\nOnly one Seurat object exists')
  stop()
}


print(SeuratOBJ.combined)
#lapply(seurat_lst, function(x) {colnames(x[[]])})
#lapply(seurat_lst, function(x) {head(x, n=3)})
#lapply(seurat_lst, function(x) max(x[["RNA"]]$counts))
#AverageExpression(SeuratOBJ, group.by = "orig.ident", features = 'ATP6AP1')

message('\nSeurats merge completed!', split(table(SeuratOBJ.combined$orig.ident), ','))

## Save merged Seurat objects

saveRDS(SeuratOBJ.combined, file = here(processedDir_out, paste0(s_sample, '.rds')))

message('Seurat merged saved!')   



## slurm script reproducibility

# slurmjobs::job_loop(
#   loops = list(type_mtx = c("data_counts", "norm_counts")),
#   name = "05_merge_seurats_job_loop_reanalyze",
#   cores = 2,
#   create_shell = TRUE
# )







############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()

