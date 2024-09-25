########################################################################
## Merge Seurat objects contained in an specific directory
##
## Authors. CSC
## Date. Feb 22, 2024
## Last.md: August 2024
##
## Input: Seurat RDS Object generated with 01_preprocessing_GEX_ATAC.R
## Output:  New Seurat merged ready to integrate with CCA/Harmony/etc
##
## NOTES: 20G free mem recommended for 20K cells
## For slurm env: runsrun --x11 --pty --partition=interactive bash
########################################################################

library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(here)

here::here()

## processed_data and plots directories 
processedDir <- here("processed-data", "02_merge_seurats")
if (!dir.exists(processedDir)) { dir.create(processedDir) }
#if (!dir.exists(here("plots", "02_merge_seurats"))) { dir.create(here("plots", "02_merge_seurats")) }
## Directory to save variable features
#dir <- here('processed-data", "02_merge_seurats", "csv_files')
#if (!dir.exists(dir)) dir.create(dir)


########################    Initials ########################  

## Counts for GEX assay available in Seurat: raw and normalized data
##        count_mtx_type_label <- 'data_counts'      
##        count_mtx_type_label <- 'norm_counts' 

args = commandArgs(trailingOnly=TRUE)
## read count mtx type:
##      all_type_mtx=(data_counts norm_counts)

count_mtx_type_label <- args[2]
message('Merging ', count_mtx_type_label, ' assays')
## for testing:  count_mtx_type_label <- "norm_counts" 

## compose rds base name for merge data
if (count_mtx_type_label=='data_counts') { 
  s_sample <- 'seurat.combined.data_counts' 
  rna_layer <- 'counts'
} else if (count_mtx_type_label=='norm_counts') { 
  s_sample <- 'seurat.combined.norm_counts' 
  rna_layer <- 'data'
} else {
  stop()
}

message('\nSuffix name for merged Seurat is `', s_sample, "`")
message('\nMerging Seurat(s) using the ', count_mtx_type_label, ' assays ')

## read directory with Seurat objects
path_directory_in <- here("processed-data", "01_preprocessing_QC")
all_rds <- paste0(path_directory_in, '/', list.files(path_directory_in, pattern="*_Hb_KDM.rds")) #, recursive = TRUE
all_rds

if (length(all_rds) < 1) { stop('\nOnly one Seurat available.') }


# ## function to save variable features before correction
# save_VFeatures <- function(sobj, f_name) {
# 
#   # Identify most highly variable genes
#   VF <- c(10,20,50,100)
#   for (x in VF) {
#     cvs_name <- ''
#     top <- head(VariableFeatures(sobj), x)
#     cvs_name <- paste0(s_sample, '_', f_name, as.character(x), '_VF.csv')
#     #print(cvs_name)
#     write.csv(top, file.path(dir, cvs_name), row.names=FALSE)
#   }
# 
# }
# 
# ## plot reductions calculated: pca, umpa, CCA and Harmony
# plot_clust <- function(sobj, f_name, reduct, ga2) {
# 
#     # integrate the samples and clusters
#     p1 <- DimPlot(sobj,
#                   reduction = reduct, group.by = c("orig.ident", ga2))
#     png_file <- paste0(f_name,'_dimplot.png')
#     png_name <- here('plots/02_merge_seurats', png_file)
#     ggsave(p1, filename = png_name, height = 5, width = 10)
# 
#     # visualize the two conditions side-by-side
#     p1 <- DimPlot(sobj,
#                   reduction = reduct, split.by = "orig.ident")
#     png_file <- paste0(f_name,'_dimplot_splitted.png')
#     png_name <- here('plots/02_merge_seurats', png_file)
#     ggsave(p1, filename = png_name, height = 5, width = 10)
# 
# }



### Prepare list of Seurat(s) to merge

seurat_lst <- list()
seurat_name_lst <- list()

message('\nPreparing Seurat objects to merge ...')



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

#seurat_lst
#seurat_name_lst
# # testing
# Cells(SeuratOBJ)[1:10]
# Features(SeuratOBJ)
# nrow(SeuratOBJ)

# some validations, 'data' slot should exists if NormalizedData was ran
# Layers(SeuratOBJ)
# max(SeuratOBJ[["RNA"]]$counts) #1255
# max(SeuratOBJ[["RNA"]]$data) #8.2943


message('\nYou merged ', length(seurat_lst), ' Seurat objects: ', sapply(seurat_name_lst, function(i) paste0(i, ' ')) )



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

saveRDS(SeuratOBJ.combined, file = here(processedDir, paste0(s_sample, '.rds')))

message('Seurat merged saved!')   



## slurm script reproducibility

# slurmjobs::job_loop(
#   loops = list(type_mtx = c("data_counts", "norm_counts")),
#   name = "01_merge_seurats_job_loop",
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

