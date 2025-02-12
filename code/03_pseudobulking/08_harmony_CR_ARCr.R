########################################################################
## Process Harmony on both GEX and ATAC 
##  
## INPUT:
##      Seurat CRr data set QCed
## OUPUT:
##      Seurat with harmonized data on both GEX and ATAC
## NOTE:
##      For +60k cells request 60G free-mem
##
## Authors. CSC 
## Date. Feb 2025
########################################################################

library("Seurat")
library("Signac") 
library("here")
#library("tidyr")
#library("stringr")

here::here()

## Check/create directories

rdsDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze")
cvsDir <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze", "cvs_files_markers")
plotsDir <- here("plots", "03_pseudobulking", "cellrangerARC_reanalyze")

if (!dir.exists(rdsDir)) {dir.create(rdsDir)}
if (!dir.exists(cvsDir)) {dir.create(cvsDir)}
if (!dir.exists(plotsDir)) {dir.create(plotsDir)}


#############################           Initials        ################################

## Set count-mtx type and integration model (CCA or Harmony)

count_mtx_type <- 'norm_counts' 
Seurat_reduction <- 'Harmony' 
# minCells <- 1

if (count_mtx_type=='data_counts') { Seurat_base_name <- 'seurat.data_counts' } else { Seurat_base_name <- 'seurat.norm_counts'}
Seurat_base_name <- here(rdsDir, paste0(Seurat_base_name, "_Harmony_ARCr_QCed.rds"))

## Validate seurat exists
if (length(list.files(rdsDir, pattern = Seurat_base_name)==0)) { message("Seurat ", basename(Seurat_base_name), " not found!"); stop() }

message("Starting cell-type identification for ", basename(Seurat_base_name))

## Load Seurat Integrated with cluster information
SeuratOBJ <- readRDS(Seurat_base_name)
## verification
# table(SeuratOBJ$orig.ident)
## Exploration
total_cells <- sum(table(SeuratOBJ$orig.ident))
# [1] 55702
Reductions(SeuratOBJ)
# [1] "pca"                "umap.unintegrated"  "integrated.cca"    
# [4] "umap"               "integrated.harmony"

## function to plot PCA and UMAP before Harmony 

f_plot_clust <- function(sobj, f_name, reduct, ga2) {
  
  # integrate the samples and clusters
  p1 <- DimPlot(sobj, 
                reduction = reduct, group.by = c("orig.ident", ga2))
  png_name <- here(plotsDir, paste0(f_name,'_dimplot.png'))  
  ggsave(p1, filename = png_name, height = 5, width = 10)
  
  # visualize the two conditions side-by-side
  p1 <- DimPlot(sobj, 
                reduction = reduct, split.by = "orig.ident")
  png_name <- here(plotsDir, paste0(f_name,'_dimplot_splitted.png'))  
  ggsave(p1, filename = png_name, height = 5, width = 10)
  
}

