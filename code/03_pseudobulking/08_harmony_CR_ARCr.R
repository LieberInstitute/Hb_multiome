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
# total_cells <- sum(table(SeuratOBJ$orig.ident))
# [1] 55702
# Reductions(SeuratOBJ)
# [1] "pca"                "umap.unintegrated"  "integrated.cca"    
# [4] "umap"               "integrated.harmony"

prefix_name <- "seurat.combined.norm_counts_rna"

p1 <- ElbowPlot(SeuratOBJ, ndims = 30, reduction = "pca") + ggtitle("Elbow on RNA") +
  geom_vline(xintercept = 20, color="red", linetype="dashed")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_elbow.png')), height = 5, width = 7)

## Re run PCA as outliers cells has been removed from the data-set (~12% = 7k)

SeuratOBJ.1 <- RunPCA(SeuratOBJ)

## some plots for comparison purposes 

p1 <- ElbowPlot(SeuratOBJ.1, ndims = 30, reduction = "pca") + ggtitle("Elbow on RNA (QCed)") +
  geom_vline(xintercept = 20, color="red", linetype="dashed") 
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_elbow_pca_QCed.png')), height = 5, width = 7)

## Re run UMAP on QCed dataset

SeuratOBJ.1 <- RunUMAP(SeuratOBJ.1, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

p1 <- DimPlot(SeuratOBJ, reduction = 'umap.unintegrated', group.by = "orig.ident") + ggtitle("UMAP on RNA")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_pca_umap.png')), height = 5, width = 7)

p1 <- DimPlot(SeuratOBJ.1, reduction = 'umap.unintegrated', group.by = "orig.ident") + ggtitle("UMAP on RNA (QCed)")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_pca_umap_QCed.png')), height = 5, width = 7)


 
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
             early_stop = T)
## rewrite harmony assay
Reductions(SeuratOBJ.1)

## Re-join layers after RNA integration
# Assays(SeuratOBJ.1)
# [1] "RNA"  "ATAC"
# SeuratOBJ.1[["RNA"]] <- JoinLayers(SeuratOBJ.1[["RNA"]])

message("Finishing Seurat-Harmony Integration on RNA - ", Sys.time())

## more plots to compare batch correction on dataset before and after remove outliers (for comparison purposes )

p1 <- DimPlot(SeuratOBJ, reduction = 'integrated.harmony', group.by = "orig.ident") + ggtitle("Harmony on RNA")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_integrated_harmony.png')), height = 5, width = 7)

p1 <- DimPlot(SeuratOBJ.1, reduction = 'integrated.harmony', group.by = "orig.ident") + ggtitle("Harmony on RNA QCed")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_integrated_harmony_QCed.png')), height = 5, width = 7)

## Run UMAP on harmonized RNA data

SeuratOBJ.1 <- RunUMAP(SeuratOBJ.1, dims = 1:30, reduction = "integrated.harmony", reduction.name = "umap.integrated")

p1 <- DimPlot(SeuratOBJ, reduction = 'umap.integrated', group.by = "orig.ident") + ggtitle("UMAP on Harmony RNA)")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_umap_integrated_harmony.png')), height = 5, width = 7)

p1 <- DimPlot(SeuratOBJ.1, reduction = 'umap.integrated', group.by = "orig.ident") + ggtitle("UMAP on Harmony RNA QCed")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_umap_integrated_harmony_QCed.png')), height = 5, width = 7)

## not longer need it
rm("SeuratOBJ")

## Run TSNE

# Reductions(SeuratOBJ.1)
SeuratOBJ.1 <- RunTSNE(SeuratOBJ.1, dims = 1:30, reduction = "integrated.harmony", reduction.name = "tsne.integrated")

p1 <- DimPlot(SeuratOBJ.1, reduction = 'tsne.integrated', group.by = "orig.ident") + ggtitle("TSNE on Harmony RNA QCed")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_TSNE_integrated_harmony_QCed.png')), height = 5, width = 7)



## Re-run TF-IDF on ATAC-QCed data

message("Run TF-IDF on ATAC")

DefaultAssay(SeuratOBJ.1) <- "ATAC"

# Run term frequency inverse document frequency (TF-IDF) normalization on a matrix.
SeuratOBJ.1 <- RunTFIDF(SeuratOBJ.1,
                        method = 1,  # computes log(𝑇𝐹×𝐼𝐷𝐹).
                        scale.factor = 10000)
SeuratOBJ.1 <- FindTopFeatures(SeuratOBJ.1,
                               min.cutoff = 'q5', # 95% most common features coverage as VariableFeatures
                               verbose = TRUE)
SeuratOBJ.1 <- RunSVD(SeuratOBJ.1)

# Reductions(SeuratOBJ.1)
# [1] "pca"                "umap.unintegrated"  "integrated.cca"    
# [4] "umap"               "integrated.harmony" "umap.integrated"   
# [7] "tsne.integrated"    "lsi"

## some plots for comparison purposes 

prefix_name <- "seurat.combined.norm_counts_atac"

p1 <- ElbowPlot(SeuratOBJ.1, ndims = 50, reduction = "lsi")  + ggtitle("Elbow on ATAC (QCed)") +
  geom_vline(xintercept = 20, color="red", linetype="dashed")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_elbow_after_QCed.png')), height = 5, width = 7)

SeuratOBJ.1 <- RunUMAP(SeuratOBJ.1, dims = 2:20, reduction = "lsi", reduction.name = "umap.lsi.unintegrated")

p1 <- DimPlot(SeuratOBJ.1, reduction = 'umap.lsi.unintegrated', group.by = "orig.ident") + ggtitle("UMAP on ATAC (QCed)")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_umap_QCed.png')), height = 5, width = 7)


message("Running Seurat-Harmony Integration on ATAC - ", Sys.time())

## split the RNA measurements into two layers one for each sample

tryCatch( { SeuratOBJ.1[["ATAC"]] <- split(SeuratOBJ.1[["ATAC"]], f = SeuratOBJ.1$orig.ident)
  return(s) }, error = function(e) { print("layers are already split") } )

## correct data on lsi by sampleID
# colnames(SeuratOBJ.1@meta.data)
SeuratOBJ.1 <- SeuratOBJ.1 |> 
  RunHarmony(group.by.vars = "orig.ident", 
             reduction.save = "integrated.lsi.harmony",
             assay.use = "ATAC",
             reduction.use= 'lsi',
             plot_convergence = TRUE,
             #max.iter = 10,  # To avoid warning() message: Quick-TRANSfer stage steps exceeded maximum (= 2785100)  
             #                  left empty as it could need more or less lters to complete
             early_stop = T,
             project.dim = F) # I think project.dim produce a bug at some point 
# `ProjectDim` arg in the harmony source code is trying to project the batch corrected embedding back to the original feature loading, 
#  set as `F` to allow the correction  

# Reductions(SeuratOBJ.1)
# [1] "pca"                    "umap.unintegrated"      "integrated.cca"        
# [4] "umap"                   "integrated.harmony"     "umap.integrated"       
# [7] "tsne.integrated"        "lsi"                    "umap.lsi.unintegrated" 
# [10] "integrated.lsi.harmony"

## Re-join layers after RNA integration
# Assays(SeuratOBJ.1)
# [1] "RNA"  "ATAC"
SeuratOBJ.1[["ATAC"]] <- JoinLayers(SeuratOBJ.1[["ATAC"]])


## more plots after correction on lsi for comparison purposes 

message("Finishing Seurat-Harmony Integration on ATAC - ", Sys.time())


## more plots after correction for comparison purposes 

p1 <- DimPlot(SeuratOBJ.1, reduction = 'integrated.lsi.harmony', group.by = "orig.ident") + ggtitle("Harmony on ATAC QCed")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_integrated_harmony_QCed.png')), height = 5, width = 7) 

## Run UMAP on harmonized RNA data
SeuratOBJ.1 <- RunUMAP(SeuratOBJ.1, dims = 2:30, reduction = "integrated.lsi.harmony", reduction.name = "umap.lsi.integrated")

p1 <- DimPlot(SeuratOBJ.1, reduction = 'umap.lsi.integrated', group.by = "orig.ident") + ggtitle("UMAP on Harmony ATAC QCed")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_umap_integrated_harmony_QCed.png')), height = 5, width = 7)

## Run TSNE
SeuratOBJ.1 <- RunTSNE(SeuratOBJ.1, dims = 2:30, reduction = "integrated.lsi.harmony", reduction.name = "tsne.lsi.integrated")

p1 <- DimPlot(SeuratOBJ.1, reduction = 'tsne.lsi.integrated', group.by = "orig.ident") + ggtitle("TSNE on Harmony ATAC QCed")
ggsave(p1, filename = here(plotsDir, paste0(prefix_name, '_TSNE_integrated_harmony_QCed.png')), height = 5, width = 7)


## Save Seurat with harmonized rna and atac data 

Seurat_base_name <- here(rdsDir, "seurat.norm_counts_ARCr_harmony_atac_rna_QCed.rds")
saveRDS(rdsDir, Seurat_base_name)

message("Saved Seurat corrected!")
