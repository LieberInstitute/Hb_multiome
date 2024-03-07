########################################################################
## Merge Seurat objects to prepare for batch effect correction
##
## Authors. CSC
## Date. Feb 22, 2024
## Last.Adaptation: xxx
##
## Input: Seurat RDS Object generated with 01_preprocessing_GEX_ATAC.R
## Output:  New Seurat combined object
##          Basic plots for reference after correction    
##
## NOTES:
## For slurm env: runsrun --x11 --pty --partition=interactive bash
########################################################################

library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
#options(Seurat.object.assay.version = 'v5')    # To use new Seurat v5: Please run: options(Seurat.object.assay.version = 'v5')
#library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
library(harmony)
#library(EnsDb.Hsapiens.v86)
#library(BSgenome.Hsapiens.UCSC.hg38)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(ggplot2)
library(patchwork)
#library(SeuratDisk)                             
library(here)

here::here()

if (!packageVersion("Seurat")=='4.9.9.9060') {
    stop
    message('This pipeline was implemented with Seurat v5 and Signac v1.11+ ')
    message('You need the laterst Seurat v5 (‘4.9.9.9060’)')
    message('Current available repository on: https://satijalab.org/seurat/articles/install.html  ') }

# Check if processed_data directory exists, if not create it
if (!dir.exists(here("processed-data/04_merge_seurats/"))) {
    dir.create(here("processed-data/04_merge_seurats/"))
}
# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/04_merge_seurats/"))) {
    dir.create(here("plots/04_merge_seurats/"))
}

source(here("code/functions_custom", "remote_plot_functions.R"))    # Call to plot GEX assay
#source(here("code/functions_custom", "remote_filtering_functions.R"))   # Call functions to subset the Seurat object

########################    Initials ########################  

## select the count-mtx to merge (raw or normalized data)
count_mtx_type <- 'data_counts'      
#count_mtx_type <- 'norm_counts' 

## load pre-existing Seurat
get_seurat <- function(name) {

        #Ex. /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/04_merge_seurats/seurat.combined.data_counts_PCA.rds
    sobj <- readRDS(name)
    # verification of the integration
    print(table(sobj$orig.ident))
    # S1_Hb_KDM S2_Hb_KDM
    # 8178      9816
    print(head(sobj, n=2))
    return(sobj)
    
}

# get QC violin plots
plot_violinQC <- function(seuratOBJ, sfeature, stitle) {
    p1 <- VlnPlot(object = seuratOBJ, features = sfeature, 
                  group.by = 'orig.ident', pt.size = 0) & geom_boxplot() &
        theme(legend.position = 'none',
              axis.text.x = element_text(angle=0, hjust=1, size=8),  #10
              axis.text.y = element_text(size=8), 
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) #&
    #labs(title = "", x = 'Samples', y ="")
    ggtitle(stitle)
    return(p1)
}

## plot reductions calculated: pca, umpa, CCA and Harmony
plot_clust <- function(sobj, f_name, reduct, ga2) {
    
    # integrate the samples and clusters
    p1 <- DimPlot(sobj, 
                  reduction = reduct, group.by = c("orig.ident", ga2))
    png_file <- paste0(f_name,'_dimplot.png')
    png_name <- here('plots/04_merge_seurats', png_file)  
    ggsave(p1, filename = png_name, height = 5, width = 10)
    
    # visualize the two conditions side-by-side
    p1 <- DimPlot(sobj, 
                  reduction = reduct, split.by = "orig.ident")
    png_file <- paste0(f_name,'_dimplot_splitted.png')
    png_name <- here('plots/04_merge_seurats', png_file)  
    ggsave(p1, filename = png_name, height = 5, width = 10)
    
}


# Directory to save variable features 
dir <- file.path(here('processed-data/04_merge_seurats/csv_files/')) 
if (!dir.exists(dir)) dir.create(dir)

## save cvs file with variable features in each reduction 
save_VFeatures <- function(sobj, f_name) {
    
    # Identify most highly variable genes
    VF <- c(10,20,50,100)
    for (x in VF) {
        cvs_name <- ''
        top <- head(VariableFeatures(sobj), x)
        cvs_name <- paste0(s_sample, '_', f_name, as.character(x), '_VF.csv')
        #print(cvs_name)
        write.csv(top, file.path(dir, cvs_name), row.names=FALSE)
    }

}





# load pre-existing seurat objects
SeuratOBJ <- get_seurat(here('processed-data/01_preprocessing_QC', 'S1_Hb_KDM.rds'))

SeuratOBJ2 <- get_seurat(here('processed-data/01_preprocessing_QC', 'S2_Hb_KDM.rds'))

# Merge Seurat objects according with the `count_mtx_type`
message('Starting merging seurats...')

if (count_mtx_type=='data_counts') {
    
    SeuratOBJ.combined <- merge(SeuratOBJ, y = SeuratOBJ2,
                                add.cell.ids = c('S1_Hb', 'S2_Hb'),
                                project = "Habenula")
    LayerData(SeuratOBJ.combined)[1:10, 1:15]
    s_sample <- 'seurat.combined.data_counts'
        
} else {
    
    SeuratOBJ <- NormalizeData(SeuratOBJ)
    SeuratOBJ2 <- NormalizeData(SeuratOBJ2)
    SeuratOBJ.combined <- merge(SeuratOBJ, y = SeuratOBJ2,
                                  add.cell.ids = c('S1_Hb', 'S2_Hb'),
                                  project = "Habenula", 
                                  merge.data = TRUE)     #  merge the normalized data matrices as well as the raw count matrices
    LayerData(SeuratOBJ.combined)[1:10, 1:15]
    s_sample <- 'seurat.combined.norm_counts'
}


#pbmc.big <- merge(pbmc3k, y = c(pbmc4k, pbmc8k), add.cell.ids = c("3K", "4K", "8K"), project = "PBMC15K")
message('Merge completed!')

# verification of the integration
table(SeuratOBJ.combined$orig.ident)
head(colnames(SeuratOBJ.combined))
tail(colnames(SeuratOBJ.combined))
unique(sapply(X = strsplit(colnames(SeuratOBJ.combined), split = "_"), FUN = "[", 1))

# > head(colnames(SeuratOBJ.combined))
# [1] "S1_Hb_AAACAGCCAAATTGCT-1" "S1_Hb_AAACAGCCAGCTTACA-1"
# [3] "S1_Hb_AAACAGCCAGGTCCTG-1" "S1_Hb_AAACAGCCAGTTAAAG-1"
# [5] "S1_Hb_AAACAGCCATAGACTT-1" "S1_Hb_AAACATGCATCCCTCA-1"
# > tail(colnames(SeuratOBJ.combined))
# [1] "S2_Hb_TTTGTGTTCGATTATG-1" "S2_Hb_TTTGTGTTCTACCTAT-1"
# [3] "S2_Hb_TTTGTTGGTAGCAGCT-1" "S2_Hb_TTTGTTGGTCGCGCAA-1"
# [5] "S2_Hb_TTTGTTGGTTTAGTCC-1" "S2_Hb_TTTGTTGGTTTGCGCC-1"

# Save RDS Object
rds_name <- here('processed-data/04_merge_seurats', paste0(s_sample, '.rds'))
# .../seurat.combined.data_counts.rds"
saveRDS(SeuratOBJ.combined, file = rds_name)
message('Seurat combined saved in ', rds_name)   


# Plot Genes and UMIs by density per cell 

# Violin plot with UMIs, Genes, ^MT and RIBO levels
p1 <- plot_violinQC(SeuratOBJ.combined, c("nCount_RNA", "nFeature_RNA", "percent.mt"), s_sample) 
#p1 <- VlnPlot(SeuratOBJ.combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "orig.ident") 
png_file <- paste0(s_sample, '_Vplots_GEX.png')
png_name <- here('plots/04_merge_seurats', png_file)  
ggsave(p1, filename = png_name, height = 5, width = 7)

# Plot Genes and UMIs by density per cell 
df_genes_per_cell <- as.data.frame(SeuratOBJ.combined[[]])
p1 <- df_genes_per_cell %>%
    ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5)) +
    #geom_vline(xintercept = 300) +
    ylab("Log10(UMIs)") +
    xlab("Gene-counts") +
    ggtitle("Genes density by cell") 

png_file <- paste0(s_sample, '_Genes_Density.png')
png_name <- here('plots/04_merge_seurats', png_file)  
ggsave(p1, filename = png_name, height = 4, width = 5)
message('UMI/Counts by MT plot saved!')  

# Gene distribution by cell
p1 <- df_genes_per_cell %>%
    ggplot(aes(x=orig.ident, y=(nFeature_RNA), fill=orig.ident)) +
    geom_boxplot(alpha = 0.7) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.text.x = element_text(vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5)) +
    ylab("Log10(nFeature_RNA)") +
    xlab("") +
    ggtitle("Genes distribution by cell")

png_file <- paste0(s_sample, '_Genes_Distribution.png')
png_name <- here('plots/04_merge_seurats', png_file)  
ggsave(p1, filename = png_name, height = 4, width = 5)
message('UMI/Counts by MT plot saved!')  

# Correlation btw Genes/UMIs 
p1 <- df_genes_per_cell %>%
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt, group.by = 'orig.ident')) + # MTRatio
    #    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=MTRatio)) + # MTRatio
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    #geom_vline(xintercept = 200, linetype=2) +
    #geom_hline(yintercept = 200, linetype=2) +
    facet_wrap(~orig.ident) +
    ylab("log10(nFeature_RNA)") +
    xlab("log10(UMIs)") +
    ggtitle('UMIs per Genes by MT levels')

png_file <- paste0(s_sample, '_UMIS_per_MT.png')
png_name <- here('plots/04_merge_seurats', png_file)  
ggsave(p1, filename = png_name, height = 4, width = 5)
message('UMI/Counts by MT plot saved!')  



############ Process PCAs

# NOTE timoast comment: https://github.com/satijalab/seurat/issues/3505
# I'd suggest doing QC and filtering cells on each object before running the integration.
# Running NormalizeData on the integrated assay will overwrite the integration results.

### Start from here / load pre-existing combined Seurat objects
if (count_mtx_type=='data_counts') { s_sample <- 'seurat.combined.data_counts' } else { s_sample <- 'seurat.combined.norm_counts' }
rds_name <- here('processed-data/04_merge_seurats', paste0(s_sample, '.rds'))
SeuratOBJ.combined <- get_seurat(rds_name)

# verification of the integration
table(SeuratOBJ.combined$orig.ident)
# S1_Hb_KDM S2_Hb_KDM 
# 8178      9816 
head(colnames(SeuratOBJ.combined))
tail(colnames(SeuratOBJ.combined))
colnames(SeuratOBJ.combined@meta.data)


######### Perform analysis without integration

# split the RNA measurements into two layers one for each sample
SeuratOBJ.combined[["RNA"]] <- split(SeuratOBJ.combined[["RNA"]], f = SeuratOBJ.combined$orig.ident)
# Warning: Assay RNA changing from Assay to Assay5
# Warning message:
#   Input is a v3 assay and `split()` only works for v5 assays; converting

# Run standard analysis: Calculate PCA cell embeddings

# LogNormalize the count data present in the assay
SeuratOBJ.combined <- NormalizeData(SeuratOBJ.combined)

# Identifies features that are outliers on a 'mean variability plot'.
# vst method (default): First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess).
SeuratOBJ.combined <- FindVariableFeatures(SeuratOBJ.combined,
                                  selection.method = "vst") # First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression...

save_VFeatures(SeuratOBJ.combined, 'pca')


# Global-scaling “LogNormalize” method that normalizes the GEX measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
all.genes <- rownames(SeuratOBJ.combined)
SeuratOBJ.combined <- ScaleData(SeuratOBJ.combined, features = all.genes)

# Run a PCA dimensionality reduction
SeuratOBJ.combined <- RunPCA(SeuratOBJ.combined)
SeuratOBJ.combined@reductions

p1 <- ElbowPlot(SeuratOBJ.combined)
png_file <- paste0(s_sample, '_PCAelbow.png')
png_name <- here('plots/04_merge_seurats', png_file)  
ggsave(p1, filename = png_name, height = 4, width = 5)


#### Before data correction some visualizations for reference
# Compute nearest neighbor graph + SNN
SeuratOBJ.combined <- FindNeighbors(SeuratOBJ.combined, dims = 1:30, reduction = "pca")
# Compute Louvain
SeuratOBJ.combined <- FindClusters(SeuratOBJ.combined, 
                          resolution = 1, cluster.name = "unintegrated_clusters")
#unique(SeuratOBJ.combined$seurat_clusters)
#unique(SeuratOBJ.combined$unintegrated_clusters)
#str(SeuratOBJ.combined)
SeuratOBJ.combined <- RunUMAP(SeuratOBJ.combined, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
head(SeuratOBJ.combined, n=2)
#SeuratOBJ.combined@reductions

# create and save UMAP-PCA plots grouped by sample and clusters and splitted side-by-side
plot_clust(SeuratOBJ.combined, paste0(s_sample, '_umap.unintegrated_'), 'umap.unintegrated', 'seurat_clusters')

# visualize more variable features in a heatmap
p1 <- DimHeatmap(SeuratOBJ.combined, reduction = 'pca', nfeatures = 30)
png_file <- paste0(s_sample, '_umap.unintegrated_pca_heatmap.png')
png_name <- here('plots/04_merge_seurats', png_file)
ggsave(p1, filename = png_name) #, height = 5, width = 5

# Save RDS Object
rds_name <- here('processed-data/04_merge_seurats', paste0(s_sample, '_PCA.rds'))
# .../seurat.combined.data_counts_PCA.rds
saveRDS(SeuratOBJ.combined, file = rds_name)
message('Seurat combined saved in ', rds_name)   



###################################################################### 
#####        Integration methods available for Seurat v5         ##### 
###################################################################### 
#############            CCAIntegration                  #############
## CCAIntegration: https://satijalab.org/seurat/reference/ccaintegration
## Example: https://satijalab.org/seurat/articles/integration_introduction.html
## “Seurat CCA” has the assumption that biologically more similar cells from different batches have a higher mathematical similarity (i.e. the dot product), and similarly, MNN assume similar cells from different batches have smaller Euclidean distance defined in the algorithm.

## Needs PCA
## The Seurat v5 integration procedure aims to return a single dimensional reduction that captures the shared sources of variance across multiple layers

### Start from here / load pre-existing combined Seurat objects
if (count_mtx_type=='data_counts') { s_sample <- 'seurat.combined.data_counts' } else { s_sample <- 'seurat.combined.norm_counts' }
rds_name <- here('processed-data/04_merge_seurats', paste0(s_sample, '_PCA.rds'))   # eurat.combined.raw_PCA.rds
#SeuratOBJ.combined <- get_seurat(rds_name)
#table(SeuratOBJ.combined$`orig.ident`)
# S1_Hb_KDM S2_Hb_KDM 
# 8178      9816 

# split the RNA measurements into two layers one for each sample, because you need to have a v5 assay (a bug ??)
#SeuratOBJ.combined[["RNA"]] <- split(SeuratOBJ.combined[["RNA"]], f = SeuratOBJ.combined$orig.ident)   #It is converted to v5 assay

message("Running Seurat-CCA Integration - ", Sys.time())
SeuratOBJ.combined <- IntegrateLayers(object = SeuratOBJ.combined, # Default dims: 1:30
                             method = CCAIntegration, 
                             orig.reduction = "pca", 
                             new.reduction = "integrated.cca",
                             verbose = FALSE)

# We can also specify parameters such as `k.anchor` to increase the strength of integration, add:  k.anchor = 20
# > SeuratOBJ.combined 
# An object of class Seurat 
# 36601 features across 17994 samples within 1 assay 
# Active assay: RNA (36601 features, 2000 variable features)
# 5 layers present: counts.S1_Hb_KDM, counts.S2_Hb_KDM, data.S1_Hb_KDM, data.S2_Hb_KDM, scale.data
# 3 dimensional reductions calculated: pca, umap.unintegrated, integrated.cca

colnames(head(SeuratOBJ.combined))

# re-join layers after integration
SeuratOBJ.combined[["RNA"]] <- JoinLayers(SeuratOBJ.combined[["RNA"]])

message("Finishing Seurat-CCA Integration - ", Sys.time())

# Cluster based in the new reduction

save_VFeatures(SeuratOBJ.combined, 'integrated.cca')

SeuratOBJ.combined <- FindNeighbors(SeuratOBJ.combined, reduction = "integrated.cca", dims = 1:30)
SeuratOBJ.combined <- FindClusters(SeuratOBJ.combined, resolution = 1)  # shared nearest neighbor (SNN)
# SeuratOBJ.combined2 <- FindClusters(SeuratOBJ.combined, resolution = 1, cluster.name = 'cca_clusters')  # shared nearest neighbor (SNN)
# Note that 'seurat_clusters' will be overwritten everytime FindClusters is run

SeuratOBJ.combined@reductions$integrated.cca
# A dimensional reduction object with key integratedcca_ 
# Number of dimensions: 50 
# Number of cells: 17994 
# Projected dimensional reduction calculated:  FALSE 
# Jackstraw run: FALSE 
# Computed using assay: RNA 

table(SeuratOBJ.combined$orig.ident)
head(SeuratOBJ.combined, n=2)

# Visualization

SeuratOBJ.combined <- RunUMAP(SeuratOBJ.combined, dims = 1:30, reduction = "integrated.cca")

# create and save UMAP-integrated.cca plots grouped by sample and clusters and splited side-by-side
plot_clust(SeuratOBJ.combined, paste0(s_sample, '_umap.integrated.cca'), 'umap', 'seurat_clusters')

# visualize more variable features in a heatmap
p1 <- DimHeatmap(SeuratOBJ.combined, reduction = 'integrated.cca', nfeatures = 30)
png_file <- paste0(s_sample, '_integrated.cca_heatmap.png')
png_name <- here('plots/04_merge_seurats', png_file)
ggsave(p1, filename = png_name, height = 5, width = 10)


rds_name <- here('processed-data/04_merge_seurats', paste0(s_sample, '_PCA_CCA.rds'))
# .../seurat.combined.data_counts_PCA_CCA.rds
saveRDS(SeuratOBJ.combined, file = rds_name)
message('Seurat combined saved in ', rds_name)   





###################################################################### 
#####        Integration methods available for Seurat v5         ##### 
###################################################################### 
##########            Harmony Integration                #############
# HarmonyIntegration: https://satijalab.org/seurat/reference/harmonyintegration
# Other options available: JointPCAIntegration, RPCAIntegration


### NOTE CSC: if correction is handled in the same Seurat integrated, then Seurat Clusters are overwriting, to keep both corrections need to be handle in separate objects 

### Start from here / load pre-existing Seurat combined objects
if (count_mtx_type=='data_counts') { s_sample <- 'seurat.combined.data_counts' } else { s_sample <- 'seurat.combined.norm_counts' }
rds_name <- here('processed-data/04_merge_seurats', paste0(s_sample, '_PCA.rds'))   # eurat.combined.raw_PCA.rds
SeuratOBJ.combined <- get_seurat(rds_name)
table(SeuratOBJ.combined$`orig.ident`)

# split the RNA measurements into two layers one for each sample, because you need to have a v5 assay (a bug ??)
#SeuratOBJ.combined[["RNA"]] <- split(SeuratOBJ.combined[["RNA"]], f = SeuratOBJ.combined$orig.ident)   #It is converted to v5 assay

message("Running Seurat-Harmony Integration - ", Sys.time())

# run correction. 
# max_iter=10 and up to 10 correction steps are expected. However, early_stop=TRUE so harmony will stop after the cost plateaus.
# Returns an object with a new dimensionality reduction
SeuratOBJ.combined <- SeuratOBJ.combined %>%
    RunHarmony(group.by.vars = "orig.ident",
               reduction = "pca",
               assay.use = 'RNA',
               reduction.save = "integrated.harmony",
               plot_convergence = TRUE,
               nclust = 50,                     # Number of clusters in model. nclust=1 equivalent to simple linear regression
               max.iter = 10,                   # One round of Harmony involves one clustering and one correction step
               #max.iter.cluster = 20,          # Maximum number of rounds to run clustering at each round of Harmony
               early_stop = T,
               #dims.use = 30,
               )
# Harmony converged after 4 iterations

#SeuratOBJ.combined@reductions

# re-join layers after integration
SeuratOBJ.combined[["RNA"]] <- JoinLayers(SeuratOBJ.combined[["RNA"]])

message("Finishing Seurat-Harmony Integration - ", Sys.time())

## Cluster based in the new reduction

save_VFeatures(SeuratOBJ.combined, 'integrated.harmony')

SeuratOBJ.combined <- FindNeighbors(SeuratOBJ.combined, reduction = "integrated.harmony", dims = 1:30)
SeuratOBJ.combined <- FindClusters(SeuratOBJ.combined, resolution = 1)
# 1 singletons identified. 17 final clusters.

SeuratOBJ.combined@reductions$integrated.harmony
# A dimensional reduction object with key harmony_ 
# Number of dimensions: 50 
# Number of cells: 17994 
# Projected dimensional reduction calculated:  FALSE 
# Jackstraw run: FALSE 
# Computed using assay: RNA 

table(SeuratOBJ.combined$orig.ident)
head(SeuratOBJ.combined, n=2)

# Visualization

SeuratOBJ.combined <- RunUMAP(SeuratOBJ.combined, dims = 1:30, reduction = "integrated.harmony")


# create and save UMAP-integrated.cca plots grouped by sample and clusters and splitted side-by-side
plot_clust(SeuratOBJ.combined, paste0(s_sample, '_umap.integrated.harmony'), 'umap', 'seurat_clusters')


# visualize more variable features in a heatmap
#p1 <- DimHeatmap(SeuratOBJ.combined, reduction = 'pca')
p1 <- DimHeatmap(SeuratOBJ.combined, reduction = 'integrated.harmony', nfeatures = 30)     # Error in Loadings(object = object, projected = projected, ...)[, dim,  :
#DimHeatmap(SeuratOBJ.combined, reduction = 'integrated.harmony', dims = 1, cells = 500, balanced = TRUE)
png_file <- paste0(s_sample, '_integrated.cca_pca_heatmap.png')
png_name <- here('plots/04_merge_seurats', png_file)
ggsave(p1, filename = png_name, height = 5, width = 10)


rds_name <- here('processed-data/04_merge_seurats', paste0(s_sample, '_PCA_Harmony.rds'))
# .../seurat.combined.data_counts_PCA_Harmony.rds
saveRDS(SeuratOBJ.combined, file = rds_name)
message('Seurat combined saved in ', rds_name)   


# INTEGRATION methods for Seurat V5:  https://satijalab.org/seurat/articles/seurat5_integration (Oct 31, 2023)
# https://satijalab.org/seurat/articles/integration_introduction.html (Nov 16, 2023)





############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()


# > Sys.time()
# [1] "2024-03-07 14:21:12 EST"
# > #"2023-04-04 12:42:26 EDT"
#   > proc.time()
# user   system  elapsed 
# 1065.858   74.277 3232.752 
# > options(width = 120)
# > session_info()
# 39m CRAN (R 4.3.2)
# data.table         1.15.0     2024-01-30 [2] CRAN (R 4.3.2)
# deldir             2.0-2      2023-11-23 [2] CRAN (R 4.3.2)
# digest             0.6.34     2024-01-11 [2] CRAN (R 4.3.2)
# dotCall64          1.1-1      2023-11-28 [2] CRAN (R 4.3.2)
# dplyr            * 1.1.4      2023-11-17 [2] CRAN (R 4.3.2)
# ellipsis           0.3.2      2021-04-29 [2] CRAN (R 4.3.2)
# fansi              1.0.6      2023-12-08 [2] CRAN (R 4.3.2)
# farver             2.1.1      2022-07-06 [2] CRAN (R 4.3.2)
# fastDummies        1.7.3      2023-07-06 [2] CRAN (R 4.3.2)
# fastmap            1.1.1      2023-02-24 [2] CRAN (R 4.3.2)
# fitdistrplus       1.1-11     2023-04-25 [2] CRAN (R 4.3.2)
# forcats          * 1.0.0      2023-01-29 [2] CRAN (R 4.3.2)
# future             1.33.1     2023-12-22 [2] CRAN (R 4.3.2)
# future.apply       1.11.1     2023-12-21 [2] CRAN (R 4.3.2)
# generics           0.1.3      2022-07-05 [2] CRAN (R 4.3.2)
# ggplot2          * 3.4.4      2023-10-12 [2] CRAN (R 4.3.2)
# ggrepel            0.9.5      2024-01-10 [2] CRAN (R 4.3.2)
# ggridges           0.5.6      2024-01-23 [2] CRAN (R 4.3.2)
# globals            0.16.2     2022-11-21 [2] CRAN (R 4.3.2)
# glue               1.7.0      2024-01-09 [2] CRAN (R 4.3.2)
# goftest            1.2-3      2021-10-07 [2] CRAN (R 4.3.2)
# gridExtra          2.3        2017-09-09 [2] CRAN (R 4.3.2)
# gtable             0.3.4      2023-08-21 [2] CRAN (R 4.3.2)
# harmony          * 1.2.0      2023-11-29 [2] CRAN (R 4.3.2)
# here             * 1.0.1      2020-12-13 [2] CRAN (R 4.3.2)
# hms                1.1.3      2023-03-21 [2] CRAN (R 4.3.2)
# htmltools          0.5.7      2023-11-03 [2] CRAN (R 4.3.2)
# htmlwidgets        1.6.4      2023-12-06 [2] CRAN (R 4.3.2)
# httpuv             1.6.14     2024-01-26 [2] CRAN (R 4.3.2)
# httr               1.4.7      2023-08-15 [2] CRAN (R 4.3.2)
# ica                1.0-3      2022-07-08 [2] CRAN (R 4.3.2)
# igraph             2.0.1.9008 2024-02-09 [2] Github (igraph/rigraph@39158c6)
# irlba              2.3.5.1    2022-10-03 [2] CRAN (R 4.3.2)
# jsonlite           1.8.8      2023-12-04 [2] CRAN (R 4.3.2)
# KernSmooth         2.23-22    2023-07-10 [3] CRAN (R 4.3.2)
# labeling           0.4.3      2023-08-29 [2] CRAN (R 4.3.2)
# later              1.3.2      2023-12-06 [2] CRAN (R 4.3.2)
# lattice            0.22-5     2023-10-24 [3] CRAN (R 4.3.2)
# lazyeval           0.2.2      2019-03-15 [2] CRAN (R 4.3.2)
# leiden             0.4.3.1    2023-11-17 [2] CRAN (R 4.3.2)
# lifecycle          1.0.4      2023-11-07 [2] CRAN (R 4.3.2)
# listenv            0.9.1      2024-01-29 [2] CRAN (R 4.3.2)
# lmtest             0.9-40     2022-03-21 [2] CRAN (R 4.3.2)
# lubridate        * 1.9.3      2023-09-27 [2] CRAN (R 4.3.2)
# magrittr           2.0.3      2022-03-30 [2] CRAN (R 4.3.2)
# MASS               7.3-60.0.1 2024-01-13 [3] CRAN (R 4.3.2)
# Matrix             1.6-5      2024-01-11 [3] CRAN (R 4.3.2)
# matrixStats        1.2.0      2023-12-11 [2] CRAN (R 4.3.2)
# mime               0.12       2021-09-28 [2] CRAN (R 4.3.2)
# miniUI             0.1.1.1    2018-05-18 [2] CRAN (R 4.3.2)
# munsell            0.5.0      2018-06-12 [2] CRAN (R 4.3.2)
# nlme               3.1-164    2023-11-27 [3] CRAN (R 4.3.2)
# parallelly         1.36.0     2023-05-26 [2] CRAN (R 4.3.2)
# patchwork        * 1.2.0      2024-01-08 [2] CRAN (R 4.3.2)
# pbapply            1.7-2      2023-06-27 [2] CRAN (R 4.3.2)
# pillar             1.9.0      2023-03-22 [2] CRAN (R 4.3.2)
# pkgconfig          2.0.3      2019-09-22 [2] CRAN (R 4.3.2)
# plotly             4.10.4     2024-01-13 [2] CRAN (R 4.3.2)
# plyr               1.8.9      2023-10-02 [2] CRAN (R 4.3.2)
# png                0.1-8      2022-11-29 [2] CRAN (R 4.3.2)
# polyclip           1.10-6     2023-09-27 [2] CRAN (R 4.3.2)
# progressr          0.14.0     2023-08-10 [2] CRAN (R 4.3.2)
# promises           1.2.1      2023-08-10 [2] CRAN (R 4.3.2)
# purrr            * 1.0.2      2023-08-10 [2] CRAN (R 4.3.2)
# R6                 2.5.1      2021-08-19 [2] CRAN (R 4.3.2)
# ragg               1.2.7      2023-12-11 [2] CRAN (R 4.3.2)
# RANN               2.6.1      2019-01-08 [2] CRAN (R 4.3.2)
# RColorBrewer       1.1-3      2022-04-03 [2] CRAN (R 4.3.2)
# Rcpp             * 1.0.12     2024-01-09 [2] CRAN (R 4.3.2)
# RcppAnnoy          0.0.22     2024-01-23 [2] CRAN (R 4.3.2)
# RcppHNSW           0.6.0      2024-02-04 [2] CRAN (R 4.3.2)
# readr            * 2.1.5      2024-01-10 [2] CRAN (R 4.3.2)
# reshape2           1.4.4      2020-04-09 [2] CRAN (R 4.3.2)
# reticulate         1.35.0     2024-01-31 [2] CRAN (R 4.3.2)
# RhpcBLASctl        0.23-42    2023-02-11 [2] CRAN (R 4.3.2)
# rlang              1.1.3      2024-01-10 [2] CRAN (R 4.3.2)
# ROCR               1.0-11     2020-05-02 [2] CRAN (R 4.3.2)
# rprojroot          2.0.4      2023-11-05 [2] CRAN (R 4.3.2)
# RSpectra           0.16-1     2022-04-24 [2] CRAN (R 4.3.2)
# Rtsne              0.17       2023-12-07 [2] CRAN (R 4.3.2)
# scales             1.3.0      2023-11-28 [2] CRAN (R 4.3.2)
# scattermore        1.2        2023-06-12 [2] CRAN (R 4.3.2)
# sctransform        0.4.1      2023-10-19 [2] CRAN (R 4.3.2)
# sessioninfo      * 1.2.2      2021-12-06 [2] CRAN (R 4.3.2)
# Seurat           * 5.0.1      2023-11-17 [2] CRAN (R 4.3.2)
# SeuratObject     * 5.0.1      2023-11-17 [2] CRAN (R 4.3.2)
# shiny              1.8.0      2023-11-17 [2] CRAN (R 4.3.2)
# sp               * 2.1-3      2024-01-30 [2] CRAN (R 4.3.2)
# spam               2.10-0     2023-10-23 [2] CRAN (R 4.3.2)
# spatstat.data      3.0-4      2024-01-15 [2] CRAN (R 4.3.2)
# spatstat.explore   3.2-6      2024-02-01 [2] CRAN (R 4.3.2)
# spatstat.geom      3.2-8      2024-01-26 [2] CRAN (R 4.3.2)
# spatstat.random    3.2-2      2023-11-29 [2] CRAN (R 4.3.2)
# spatstat.sparse    3.0-3      2023-10-24 [2] CRAN (R 4.3.2)
# spatstat.utils     3.0-4      2023-10-24 [2] CRAN (R 4.3.2)
# stringi            1.8.3      2023-12-11 [2] CRAN (R 4.3.2)
# stringr          * 1.5.1      2023-11-14 [2] CRAN (R 4.3.2)
# survival           3.5-7      2023-08-14 [3] CRAN (R 4.3.2)
# systemfonts        1.0.5      2023-10-09 [2] CRAN (R 4.3.2)
# tensor             1.5        2012-05-05 [2] CRAN (R 4.3.2)
# textshaping        0.3.7      2023-10-09 [2] CRAN (R 4.3.2)
# tibble           * 3.2.1      2023-03-20 [2] CRAN (R 4.3.2)
# tidyr            * 1.3.1      2024-01-24 [2] CRAN (R 4.3.2)
# tidyselect         1.2.0      2022-10-10 [2] CRAN (R 4.3.2)
# tidyverse        * 2.0.0      2023-02-22 [2] CRAN (R 4.3.2)
# timechange         0.3.0      2024-01-18 [2] CRAN (R 4.3.2)
# tzdb               0.4.0      2023-05-12 [2] CRAN (R 4.3.2)
# utf8               1.2.4      2023-10-22 [2] CRAN (R 4.3.2)
# uwot               0.1.16     2023-06-29 [2] CRAN (R 4.3.2)
# vctrs              0.6.5      2023-12-01 [2] CRAN (R 4.3.2)
# viridisLite        0.4.2      2023-05-02 [2] CRAN (R 4.3.2)
# withr              3.0.0      2024-01-16 [2] CRAN (R 4.3.2)
# xtable             1.8-4      2019-04-21 [2] CRAN (R 4.3.2)
# zoo                1.8-12     2023-04-13 [2] CRAN (R 4.3.2)
# 
# [1] /users/csoto/R/4.3.x
# [2] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3.x/R/lib64/R/library
