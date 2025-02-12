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
library("ggplot2")
library("harmony")
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


prefix_name <- "seurat.combined.norm_counts"

p1 <- ElbowPlot(SeuratOBJ, ndims = 30, reduction = "pca")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_elbow_before_QCed.png')), height = 5, width = 7)

## Re run PCA as outliers cells has been removed from the data-set (~12% = 7k)
SeuratOBJ.1 <- RunPCA(SeuratOBJ)
p1 <- ElbowPlot(SeuratOBJ.1, ndims = 30, reduction = "pca")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_elbow_after_QCed.png')), height = 5, width = 7)


## Re run UMAP
SeuratOBJ.1 <- RunUMAP(SeuratOBJ.1, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

p1 <- DimPlot(SeuratOBJ, reduction = 'umap.unintegrated', group.by = "orig.ident")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_umap.png')), height = 5, width = 7)

p1 <- DimPlot(SeuratOBJ.1, reduction = 'umap.unintegrated', group.by = "orig.ident")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_umap_QCed.png')), height = 5, width = 7)


 
## Re-run Harmony on QCed data
## NOTE: if correction is handled in the same Seurat integrated, then Seurat Clusters are overwriting

message("Running Seurat-Harmony Integration on RNA - ", Sys.time())

# dd/mm/yyyy
set.seed(12022025)

# use k-means centroids initialization
SeuratOBJ.1 <- SeuratOBJ.1 |>
  RunHarmony(group.by.vars = "orig.ident",
             reduction = "pca",
             assay.use = "RNA",
             reduction.save = "integrated.harmony",
             plot_convergence = TRUE,
             #nclust = 50,                     # Number of clusters in model. nclust=1 equivalent to simple linear regression
             max.iter = 10,                   # One round of Harmony involves one clustering and one correction step
             #max.iter.cluster = 20,          # Maximum number of rounds to run clustering at each round of Harmony
             early_stop = T
  )
## rewrite harmony assay
Reductions(SeuratOBJ.1)

message("Finishing Seurat-Harmony Integration on RNA - ", Sys.time())

p1 <- DimPlot(SeuratOBJ, reduction = 'integrated.harmony', group.by = "orig.ident")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_umap_integrated_harmony.png')), height = 5, width = 7)

p1 <- DimPlot(SeuratOBJ.1, reduction = 'integrated.harmony', group.by = "orig.ident")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_umap_integrated_harmony_QCed.png')), height = 5, width = 7)


## Run UMAP on harmonized RNA data
SeuratOBJ.1 <- RunUMAP(SeuratOBJ.1, dims = 1:30, reduction = "integrated.harmony", reduction.name = "umap.integrated")

p1 <- DimPlot(SeuratOBJ, reduction = 'integrated.harmony', group.by = "orig.ident")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_umap_integrated_harmony.png')), height = 5, width = 7)

p1 <- DimPlot(SeuratOBJ.1, reduction = 'integrated.harmony', group.by = "orig.ident")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_umap_integrated_harmony_QCed.png')), height = 5, width = 7)


## Run TSNE
Reductions(SeuratOBJ.1)
SeuratOBJ.1 <- RunTSNE(SeuratOBJ.1, dims = 1:30, reduction = "integrated.harmony", reduction.name = "tsne.integrated")

p1 <- DimPlot(SeuratOBJ.1, reduction = 'tsne.integrated', group.by = "orig.ident")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_TSNE_integrated_harmony_QCed.png')), height = 5, width = 7)




