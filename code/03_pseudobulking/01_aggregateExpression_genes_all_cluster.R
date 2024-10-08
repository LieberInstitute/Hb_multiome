########################################################################

## Aggregate gene expression (from Seurat) for RNA assay counts and data (normalized) for CCA and Harmony reductions
## Authors. CSC/lcollado
## Last md: Aug, 2024
##
## Input:  Seurat `CCA` or `Harmony integration RDS object
## Output:  (1) Seurat pseudo bulked 
##          (2) DEG cvs file before pseudo bulk 
##          (3) DEG cvs file after pseudo bulk

## AggregateExpression passes inputs to PseudobulkExpression. The outputs will be the same assuming the input parameters are identical.
## When running AggregateExpression on an integrated assay then it should reflect any batch correction that was performed, assuming you've specified the correct assays value. AggregateExpression does not perform any batch correction itself.
## AggregateExpression is only intended to be run on the raw counts. You could call PseudobulkExpression using method="aggregate".
##
## NOTES: recommended ~80G free-mem for 50-70G cells

########################################################################

library("Seurat")                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
## required packages for aggregation
library("multtest")
library("metap")
library("tidyverse")
library("here")

here::here()

########################    Initials ######################## 

# Check/create directories 
inputDir <- here("processed-data", "02_merge_seurats", "cellrangerARC_reanalyze")
processedDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze")
plotDir <- here("plots", "03_pseudobulking", "cellrangerARC_reanalyze")
cvsDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze", "cvs_files_markers")

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }
if (!dir.exists(cvsDir)) { dir.create(cvsDir) }

## For testing: By Default `norm_counts Harmony`
## For index 2:
#   count_mtx_type <- 'data_counts'
#   count_mtx_type <- 'norm_counts' 
## For index 4:
#   integration_model <- 'CCA'
#   integration_model <- 'Harmony'

## read input arguments
args = commandArgs(trailingOnly=TRUE)
count_mtx_type <- args[2]
integration_model <- args[4]

if (length(count_mtx_type) && length(integration_model)) {
  message("Processing ", count_mtx_type, " with ", integration_model, " method.")
} else {
  message("Input arguments missed")
  message("Assay type: ", count_mtx_type, " Intergration model: ", integration_model)
  stop()
}

## function to load pre-existing Seurat
get_seurat <- function(name) { sobj <- readRDS(name); return(sobj)}

## Build Seurat object name
if (count_mtx_type=='data_counts') { 
  Seurat_base_name <- 'seurat.data_counts' 
} else { 
  Seurat_base_name <- 'seurat.norm_counts' 
}
if (integration_model=='CCA') {
  rds_name <- here(inputDir, paste0(Seurat_base_name, '_CCA.rds'))
} else {
  rds_name <- here(inputDir, paste0(Seurat_base_name, '_Harmony.rds'))
}
#rds_name

message("Processing ", Seurat_base_name, " ", integration_model)

## Load Seurat object
SeuratOBJ <- get_seurat(rds_name)

## Verify Seurats to aggregate
table(SeuratOBJ$orig.ident)
# head(colnames(SeuratOBJ))
# tail(colnames(SeuratOBJ))
# colnames(SeuratOBJ@meta.data)
SeuratOBJ@reductions
table(SeuratOBJ$seurat_clusters)

###################### Pseudo bulk expression data  ######################

## Run in an integrated Seurat
SeuratOBJ[["RNA"]] <- JoinLayers(SeuratOBJ[["RNA"]])

## count cells by clusters
unique(SeuratOBJ@meta.data$seurat_clusters)
table(Idents(SeuratOBJ))

## Find DEG in the integrated Seurat for ALL clusters (BEFORE pseudobulk)
all.markers <- FindAllMarkers(object = SeuratOBJ)
#head(all.markers, n=3)

# cvs_file <- paste0(Seurat_base_name, '_', integration_model, '_Allmarkers.csv')
cvs_file <- paste0(Seurat_base_name, '_', integration_model, '_Allmarkers.csv')
cvs_file <- here(cvsDir, cvs_file)
write.csv(all.markers, cvs_file)

message(" FindAllMarkers in batch corrected data done!")

## Save new Seurat pseudo bulk 
rds_name <- paste0(Seurat_base_name,'_', integration_model, '_All.rds')
rds_name <- here(processedDir, rds_name)
saveRDS(SeuratOBJ, file = rds_name)

message(" Saved Seurat batch corrected data clusters not filtered.")

## Apply pseudo bulk to ALL clusters and marker genes selected. To remove cluster with few cells update 01_aggregate_GeneExpression_filtered_clusters.R script 
## AggregateExpression from Seurat average data.  But can be used this internal Seurat:::PseudobulkExpression(pb.method = 'aggregate' ) to sum up counts by `categories
SeuratOBJ_Hb_all_pseudobulked <- Seurat:::AggregateExpression(SeuratOBJ, 
                                                     return.seurat = TRUE,
                                                     group.by = c("seurat_clusters", "orig.ident"))

## https://github.com/satijalab/seurat/issues/8919#issuecomment-2125129658 
# ## AggregateExpression with return.seurat=FALSE will return the summed counts

## If return.seurat = TRUE, aggregated values are placed in the 'counts' layer of the returned object. The data is then normalized by running NormalizeData on the aggregated counts. ScaleData is then run on the default assay before returning the object.

SeuratOBJ_Hb_all_pseudobulked

table(SeuratOBJ_Hb_all_pseudobulked$seurat_clusters)
# g0  g1 g10 g11 g12 g13 g14 g15 g17 g18 g19  g2 g21  g3  g4  g5  g6  g7  g9 
# 3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3   3

all.markers_p <- FindAllMarkers(object = SeuratOBJ_Hb_all_pseudobulked)
# Warning: When testing g18_S2-Hb-KDM versus all:
#     Cell group 1 has fewer than 3 cells

## check DEG found in the pseudo bulk data
if (length(all.markers_p)>0) {
  cvs_file <- paste0(Seurat_base_name, '_', integration_model, '_Allmarkers_cells_pseudobulk.csv')
  cvs_file <- here(cvsDir, cvs_file)
  write.csv(all.markers_p, cvs_file)
} else {
  message("None FindAllMarkers() found in pseudobulk data.")
}

## Save new Seurat pseudo bulk 
rds_name <- paste0(Seurat_base_name,'_', integration_model, '_All_pseudobulk.rds')
rds_name <- here(processedDir, rds_name)
saveRDS(SeuratOBJ_Hb_all_pseudobulked, file = rds_name)

message(" Pseudobulk completed!")


## slurm script reproducibility

# slurmjobs::job_loop(
#   loops = list(type_mtx = c("data_counts", "norm_counts"), integration_model = c("CCA", "Harmony")),
#   name = "01_aggregateExpression_genes_all_cluster",
#   cores = 2,
#   create_shell = TRUE
# )


############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()

