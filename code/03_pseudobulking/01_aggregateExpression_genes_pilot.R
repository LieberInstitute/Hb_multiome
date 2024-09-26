########################################################################

## Pseudo bulk with Aggregate expression (from Seurat) for RNA assay for CCA and Harmony reductions
## Authors. CSC/lcollado
## Last md: Aug, 2024
##
## Input:  Seurat integrated object with integrated samples. Ex. Habenula samples S1 and S2 
## Output:  (1) Seurat pseudobulked, 
##          (2) DEG cvs file before pseudobulk, 
##          (3) DEG cvs file after pseudobulk

## AggregateExpression passes inputs to PseudobulkExpression. The outputs will be the same assuming the input parameters are identical.
## When running AggregateExpression on an integrated assay then it should reflect any batch correction that was performed, assuming you've specified the correct assays value. AggregateExpression does not perform any batch correction itself.
## AggregateExpression is only intended to be run on the raw counts. You could call PseudobulkExpression using method="aggregate".
##
## NOTES: recommended ~20G free-mem

########################################################################

library('Seurat')                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
## required packages for aggregation
library('multtest')
library('metap')
library('tidyverse')
# library('ggplot2')
# library('patchwork')
library(here)

here::here()

########################    Initials ######################## 

# Check/create directories 
processedDir <- here("processed-data", "03_pseudobulking", "pilot_results")
if (!dir.exists(processedDir)) { dir.create(processedDir) }
plotDir <- here("plots", "03_pseudobulking", "pilot_results")
if (!dir.exists(plotDir)) { dir.create(plotDir) }
cvsDir <- here("processed-data", "03_pseudobulking", "cvs_files_markers", "pilot_results")
if (!dir.exists(cvsDir)) { dir.create(cvsDir) }
inputDir <- here("processed-data", "02_merge_seurats", "pilot_results") 

## Select the count-mtx to merge (raw or normalized data). By Default `norm_counts Harmony`
count_mtx_type <- 'data_counts'
#count_mtx_type <- 'norm_counts' 
#seurat_reduction <- 'CCA'
Seurat_reduction <- 'Harmony'

## function to load pre-existing Seurat
get_seurat <- function(name) { sobj <- readRDS(name); return(sobj)}

## Compose Seurat object name
if (count_mtx_type=='data_counts') { 
  Seurat_base_name <- 'seurat.combined.data_counts' 
} else { 
  Seurat_base_name <- 'seurat.combined.norm_counts' 
}
if (Seurat_reduction=='CCA') {
  rds_name <- here(inputDir, paste0(Seurat_base_name, '_PCA_CCA.rds'))
} else {
  rds_name <- here(inputDir, paste0(Seurat_base_name, '_PCA_Harmony.rds'))
}
rds_name
## Set minimum number of cells by cluster to process
min_cells <- 1

## Load Seurat object
SeuratOBJ <- get_seurat(rds_name)
# An object of class Seurat 
# 36601 features across 17994 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
# 3 layers present: data, counts, scale.data
# 4 dimensional reductions calculated: pca, umap.unintegrated, integrated.harmony, umap

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

#### Calculate the minimum number of cells by cluster.
##      To avoid issue when having few cells need to adjust the minimum number of cells
##      For example, if there is a cluster "15" that has 0 cells, the function will skip that cluster with a warning (that's perfect). Also if the number of cells is between min.cells.groups (default = 1-3) and 0, an error is thrown and it stops working. That is why I previously remove from the Seurat Object the cells of the clusters with 3 or less cells for each condition/sample. 

few_cells_samples <- unique(SeuratOBJ@meta.data$orig.ident)
few_cells <- vector()
few_cells_tmp <- vector()
## Remove clusters with a minimum number of cells `min_cells` in each sample/condition. Usually from 1 to 3 cell.
for (i in 1:length(few_cells_samples)) {    
  # for testing: i <- '4S_Hb_KDM'
  message("Sample: ", i, " ")
  table(SeuratOBJ@meta.data$seurat_clusters[SeuratOBJ@meta.data$orig.ident == few_cells_samples[i]])
  few_cells_tmp <- table(SeuratOBJ@meta.data$seurat_clusters[SeuratOBJ@meta.data$orig.ident == few_cells_samples[i]]) <= min_cells
  few_cells_tmp <- names(few_cells_tmp)[few_cells_tmp == "TRUE"]
  few_cells <- c(few_cells,few_cells_tmp)
}
message(" Clusters with less than ", min_cells," cell: ", length(unique(few_cells)))
print(table(few_cells))

## Subset and keep clusters with cells => `min_cells`
clusters <- sort(unique(SeuratOBJ@meta.data$seurat_clusters))
`%notin%` <- Negate(`%in%`) 
clusters_high_cell <- clusters[clusters %notin% unique(few_cells)]
SeuratOBJ <- subset(SeuratOBJ, subset = seurat_clusters %in% as.vector(clusters_high_cell))
#SeuratOBJtmp<-SeuratOBJ
## count cells by clusters
unique(SeuratOBJ@meta.data$seurat_clusters)
table(Idents(SeuratOBJ))

## Find DEG in the integrated Seurat for ALL clusters (BEFORE pseudobulk)
all.markers <- FindAllMarkers(object = SeuratOBJ)
#head(all.markers, n=3)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster   gene
# NPAS3      0  -1.836151 0.078 0.736         0       0  NPAS3

# cvs_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_Allmarkers.csv')
cvs_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_Allmarkers_min', min_cells, 'cells_pilot.csv')
cvs_file <- here(cvsDir, cvs_file)
write.csv(all.markers, cvs_file)

message(" FindAllMarkers in batch corrected data done!")

## Save new Seurat pseudo bulk 
rds_name <- paste0(Seurat_base_name,'_', Seurat_reduction, '_subset_pilot.rds')
rds_name <- here(processedDir, rds_name)
saveRDS(SeuratOBJ, file = rds_name)

message(" Saved Seurat batch corrected data with `min_cell=", min_cells,"` clusters filtered")



## Apply pseudo bulk to ALL clusters and marker genes selected
# ~/seurat.combined.data_counts_PCA_Harmony_DoHeatmap_pseudobulk.pdf
## Currently AggregateExpression can only average your data.  But you can use this internal Seurat:::PseudobulkExpression(pb.method = 'aggregate' ) to sum up counts by `categories
SeuratOBJ_Hb_all_pseudobulked <- Seurat:::AggregateExpression(SeuratOBJ, 
                                                     return.seurat = TRUE,
                                                     group.by = c("seurat_clusters", "orig.ident"))

## https://github.com/satijalab/seurat/issues/8919#issuecomment-2125129658 
# ## AggregateExpression with return.seurat=FALSE will return the summed counts

## If return.seurat = TRUE, aggregated values are placed in the 'counts' layer of the returned object. The data is then normalized by running NormalizeData on the aggregated counts. ScaleData is then run on the default assay before returning the object.

SeuratOBJ_Hb_all_pseudobulked

table(SeuratOBJ_Hb_all_pseudobulked$seurat_clusters)
# g0  g1 g10 g11 g12 g13 g14 g15  g2  g3  g4  g5  g6  g7  g8  g9 
# 2   2   2   2   2   2   2   2   2   2   2   2   2   2   2   2

all.markers_p <- FindAllMarkers(object = SeuratOBJ_Hb_all_pseudobulked)
# Warning: When testing g18_S2-Hb-KDM versus all:
#     Cell group 1 has fewer than 3 cells

## check DEG found in the pseudobulk data
if (length(all.markers_p)>0) {
  cvs_file <- paste0(Seurat_base_name, '_', Seurat_reduction, '_Allmarkers_min', min_cells, 'cells_pseudobulk.csv')
  cvs_file <- here(cvsDir, cvs_file)
  write.csv(all.markers_p, cvs_file)
} else {
  message("None FindAllMarkers() with more then ", min_cells, " found in pseudobulk data.")
}

## Save new Seurat pseudo bulk 
rds_name <- paste0(Seurat_base_name,'_', Seurat_reduction, '_pseudobulk_pilot.rds')
rds_name <- here(processedDir, rds_name)
saveRDS(SeuratOBJ_Hb_all_pseudobulked, file = rds_name)

message(" Pseudobulk done!")


############################ Move plots to 02_aggregateExpression_plots.R script. CSC ############################


# ## slurm script reproducibility
# 
# slurmjobs::job_loop(
#   loops = list(type_mtx = c("data_counts", "norm_counts"), integration_model = c("mod1", "mod2", "mod3", "mod4")),
#   name = "01_aggregateExpression_genes",
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

