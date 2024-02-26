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
## For slurm env: $srun --x11 --pty --partition=interactive bash
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
if (!dir.exists(here("processed-data/01_preprocessing_QC/"))) {
    dir.create(here("processed-data/01_preprocessing_QC/"))
}
# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/01_preprocessing_QC/"))) {
    dir.create(here("plots/01_preprocessing_QC/"))
}

source(here("code/functions_custom", "remote_plot_functions.R"))    # Call to plot GEX assay
#source(here("code/functions_custom", "remote_filtering_functions.R"))   # Call functions to subset the Seurat object

########################    Initials ########################  

## select the count-mtx to merge (raw or normalized data)
count_mtx_type <- 'raw_counts'      
#count_mtx_type <- 'normalized' 

## load pre-existing Seurat
get_seurat <- function(name) {

        #Ex. /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/01_preprocessing_QC/seurat.combined.raw_PCA.rds
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
plot_clust <- function(reduct) {
    # reduct <-'pca'
    # integrate the samples and clusters
    p1 <- DimPlot(SeuratOBJ, 
                  reduction = reduct, group.by = c("orig.ident", "seurat_clusters"))
    png_file <- paste0(s_sample, '_',reduct,'_dimplot.png')
    png_name <- here('plots/01_preprocessing_QC', png_file)  
    ggsave(p1, filename = png_name, height = 5, width = 10)
    
    # visualize the two conditions side-by-side
    p1 <- DimPlot(SeuratOBJ, 
                  reduction = reduct, split.by = "orig.ident")
    png_file <- paste0(s_sample, '_',reduct,'_dimplot_splitted.png')
    png_name <- here('plots/01_preprocessing_QC', png_file)  
    ggsave(p1, filename = png_name, height = 5, width = 10)
    
    # # visualize more variable features in a heatmap
    # p1 <- DimHeatmap(SeuratOBJ,
    #                  reduction = reduct, nfeatures = 30)
    # png_file <- paste0(s_sample, '_',reduct,'_Heatmap.png')
    # png_name <- here('plots/01_preprocessing_QC', png_file)  
    # ggsave(p1, filename = png_name, height = 5, width = 10)
    
}


# Directory to save variable features 
dir <- file.path(here('processed-data/01_preprocessing_QC/csv_files/')) 
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

if (count_mtx_type=='raw_counts') {
    
    SeuratOBJ.combined <- merge(SeuratOBJ, y = SeuratOBJ2,
                                add.cell.ids = c('S1_Hb_KDM', 'S2_Hb_KDM'),
                                project = "Habenula")
    LayerData(SeuratOBJ.combined)[1:10, 1:15]
    s_sample <- 'seurat.combined.raw'
        
} else {
    
    SeuratOBJ <- NormalizeData(SeuratOBJ)
    SeuratOBJ2 <- NormalizeData(SeuratOBJ2)
    SeuratOBJ.combined <- merge(SeuratOBJ, y = SeuratOBJ2,
                                  add.cell.ids = c('S1_Hb_KDM', 'S2_Hb_KDM'),
                                  project = "Habenula", 
                                  merge.data = TRUE)     #  merge the normalized data matrices as well as the raw count matrices
    LayerData(SeuratOBJ.combined)[1:10, 1:15]
    s_sample <- 'seurat.combined.normalized'
}


#pbmc.big <- merge(pbmc3k, y = c(pbmc4k, pbmc8k), add.cell.ids = c("3K", "4K", "8K"), project = "PBMC15K")
message('Merge completed!')

# verification of the integration
table(SeuratOBJ.combined$orig.ident)
head(colnames(SeuratOBJ.combined))
tail(colnames(SeuratOBJ.combined))
unique(sapply(X = strsplit(colnames(SeuratOBJ.combined), split = "_"), FUN = "[", 1))

# > head(colnames(SeuratOBJ.combined))
# [1] "S1_Hb_KDM_AAACAGCCAAATTGCT-1" "S1_Hb_KDM_AAACAGCCAGCTTACA-1"
# [3] "S1_Hb_KDM_AAACAGCCAGGTCCTG-1" "S1_Hb_KDM_AAACAGCCAGTTAAAG-1"
# [5] "S1_Hb_KDM_AAACAGCCATAGACTT-1" "S1_Hb_KDM_AAACATGCATCCCTCA-1"
# > tail(colnames(SeuratOBJ.combined))
# [1] "S2_Hb_KDM_TTTGTGTTCGATTATG-1" "S2_Hb_KDM_TTTGTGTTCTACCTAT-1"
# [3] "S2_Hb_KDM_TTTGTTGGTAGCAGCT-1" "S2_Hb_KDM_TTTGTTGGTCGCGCAA-1"
# [5] "S2_Hb_KDM_TTTGTTGGTTTAGTCC-1" "S2_Hb_KDM_TTTGTTGGTTTGCGCC-1"

# Save RDS Object
rds_name <- here('processed-data/01_preprocessing_QC', paste0(s_sample, '.rds'))
# [1] "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/01_preprocessing_QC/seurat.combined.raw.rds"
saveRDS(SeuratOBJ.combined, file = rds_name)
message('Seurat combined saved in ', rds_name)   


# Plot Genes and UMIs by density per cell 

# Violin plot with UMIs, Genes, ^MT and RIBO levels
p1 <- plot_violinQC(SeuratOBJ.combined, c("nCount_RNA", "nFeature_RNA", "percent.mt"), s_sample) 
#p1 <- VlnPlot(SeuratOBJ.combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "orig.ident") 
png_file <- paste0(s_sample, '_Vplots_GEX.png')
png_name <- here('plots/01_preprocessing_QC', png_file)  
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
png_name <- here('plots/01_preprocessing_QC', png_file)  
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
png_name <- here('plots/01_preprocessing_QC', png_file)  
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
png_name <- here('plots/01_preprocessing_QC', png_file)  
ggsave(p1, filename = png_name, height = 4, width = 5)
message('UMI/Counts by MT plot saved!')  



############ Process PCAs

# NOTE timoast comment: https://github.com/satijalab/seurat/issues/3505
# I'd suggest doing QC and filtering cells on each object before running the integration.
# Running NormalizeData on the integrated assay will overwrite the integration results.

### Start from here / load pre-existing seurat objects
if (count_mtx_type=='raw_counts') { s_sample <- 'seurat.combined.raw' } else { s_sample <- 'seurat.combined.normalized' }
rds_name <- here('processed-data/01_preprocessing_QC', paste0(s_sample, '.rds'))
SeuratOBJ <- get_seurat(rds_name)

######### Perform analysis without integration

# split the RNA measurements into two layers one for each sample
SeuratOBJ.combined[["RNA"]] <- split(SeuratOBJ.combined[["RNA"]], f = SeuratOBJ.combined$orig.ident)


SeuratOBJ <- SeuratOBJ.combined     # Warning: Assay RNA changing from Assay to Assay5

# Run standard analysis: Calculate PCA cell embeddings

# LogNormalize the count data present in the assay
SeuratOBJ <- NormalizeData(SeuratOBJ)

# Identifies features that are outliers on a 'mean variability plot'.
# vst method (default): First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess).
SeuratOBJ <- FindVariableFeatures(SeuratOBJ,
                                  selection.method = "vst") # First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression...

save_VFeatures(SeuratOBJ, 'pca')


# Global-scaling “LogNormalize” method that normalizes the GEX measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
all.genes <- rownames(SeuratOBJ)
SeuratOBJ <- ScaleData(SeuratOBJ, features = all.genes)

# Run a PCA dimensionality reduction
SeuratOBJ <- RunPCA(SeuratOBJ)

p1 <- ElbowPlot(SeuratOBJ)
png_file <- paste0(s_sample, '_PCAelbow.png')
png_name <- here('plots/01_preprocessing_QC', png_file)  
ggsave(p1, filename = png_name, height = 4, width = 5)

# Save RDS Object
rds_name <- here('processed-data/01_preprocessing_QC', paste0(s_sample, '_PCA.rds'))
# .../seurat.combined.raw.rds"
saveRDS(SeuratOBJ, file = rds_name)
message('Seurat combined saved in ', rds_name)   


#### Before data correction some visualizations for reference

SeuratOBJ <- FindNeighbors(SeuratOBJ, dims = 1:30, reduction = "pca")
SeuratOBJ <- FindClusters(SeuratOBJ, 
                          resolution = 2, cluster.name = "unintegrated_clusters")
SeuratOBJ <- RunUMAP(SeuratOBJ, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
head(SeuratOBJ, n=2)
#SeuratOBJ@reductions

Map(plot_clust, 'umap.unintegrated')



###################################################################### 
#####        Integration methods available for Seurat v5         ##### 
###################################################################### 
#############            CCAIntegration                  #############
## CCAIntegration: https://satijalab.org/seurat/reference/ccaintegration
## “Seurat CCA” has the assumption that biologically more similar cells from different batches have a higher mathematical similarity (i.e. the dot product), and similarly, MNN assume similar cells from different batches have smaller Euclidean distance defined in the algorithm.

## Needs PCA

### Start from here / load pre-existing Seurat objects
if (count_mtx_type=='raw_counts') { s_sample <- 'seurat.combined.raw' } else { s_sample <- 'seurat.combined.normalized' }
rds_name <- here('processed-data/01_preprocessing_QC', paste0(s_sample, '_PCA.rds'))
SeuratOBJ <- get_seurat(rds_name)
table(SeuratOBJ$`orig.ident`)
# S1_Hb_KDM S2_Hb_KDM 
# 8178      9816 

# split the RNA measurements into two layers one for each sample, because you need to have a v5 assay (a bug ??)
SeuratOBJ[["RNA"]] <- split(SeuratOBJ[["RNA"]], f = SeuratOBJ$orig.ident)   #It is converted to v5 assay

message("Running Seurat-CCA Integration - ", Sys.time())
SeuratOBJ <- IntegrateLayers(object = SeuratOBJ, # Default dims: 1:30
                             method = CCAIntegration, 
                             orig.reduction = "pca", 
                             new.reduction = "integrated.cca",
                             verbose = FALSE)

# We can also specify parameters such as `k.anchor` to increase the strength of integration, add:  k.anchor = 20

# re-join layers after integration
SeuratOBJ[["RNA"]] <- JoinLayers(SeuratOBJ[["RNA"]])

message("Finishing Seurat-CCA Integration - ", Sys.time())

## Run standard analysis: Calculate integrated.cca cell embeddings

# LogNormalize the count data present in the assay
SeuratOBJ <- NormalizeData(SeuratOBJ)

# Identifies features that are outliers on a 'mean variability plot'.
# vst method (default): First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess).
SeuratOBJ <- FindVariableFeatures(SeuratOBJ,
                                  selection.method = "vst") # First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression...

map(save_VFeatures, 'integrated.cca')
# # Scale data
# # Run a PCA dimensionality reduction

SeuratOBJ <- FindNeighbors(SeuratOBJ, reduction = "integrated.cca", dims = 1:30)
SeuratOBJ <- FindClusters(SeuratOBJ, resolution = 1)

SeuratOBJ <- RunUMAP(SeuratOBJ, dims = 1:30, reduction = "integrated.cca")

SeuratOBJ@reductions$integrated.cca
# A dimensional reduction object with key integratedcca_ 
# Number of dimensions: 50 
# Number of cells: 17994 
# Projected dimensional reduction calculated:  FALSE 
# Jackstraw run: FALSE 
# Computed using assay: RNA 

table(SeuratOBJ$orig.ident)
head(SeuratOBJ, n=2)

## Visualization
Map(plot_clust, 'integrated.cca')

rds_name <- here('processed-data/01_preprocessing_QC', paste0(s_sample, '_PCA_CCA.rds'))
# [1] "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/01_preprocessing_QC/seurat.combined.raw.rds"
saveRDS(SeuratOBJ, file = rds_name)
message('Seurat combined saved in ', rds_name)   

# # Directory to save variable features 
# dir <- file.path(here('processed-data/01_preprocessing_QC/csv_files/')) 
# if (!dir.exists(dir)) dir.create(dir)
# 
# # Identify most highly variable genes
# VF <- c(10,20,50,100)
# 
# for (x in VF) {
#     top <- head(VariableFeatures(SeuratOBJ), x)
#     write.csv(top, file.path(dir, paste0(s_sample,'_', x,'_CCAIntegrated_VF.csv')), row.names=FALSE)
# }



###################################################################### 
#####        Integration methods available for seurat v5         ##### 
###################################################################### 
##########            Harmony Integration                #############
# HarmonyIntegration: https://satijalab.org/seurat/reference/harmonyintegration
# Other options available: JointPCAIntegration, RPCAIntegration

if (count_mtx_type=='raw_counts') { s_sample <- 'seurat.combined.raw' } else { s_sample <- 'seurat.combined.normalized' }
rds_name <- here('processed-data/01_preprocessing_QC', paste0(s_sample, '_PCA_CCA.rds'))
SeuratOBJ <- get_seurat(rds_name)

# max_iter=10 and up to 10 correction steps are expected. However, early_stop=TRUE so harmony will stop after the cost plateaus.
# Returns an object with a new dimensionality reduction
message("Running Seurat-Harmony Integration - ", Sys.time())

SeuratOBJ <- SeuratOBJ %>% 
    RunHarmony(group.by.vars = "orig.ident", 
               reduction = "pca",
               assay.use = 'RNA',
               reduction.save = "harmony_r",
               plot_convergence = TRUE, 
               nclust = 50,                    # Number of clusters in model. nclust=1 equivalent to simple linear regression
               max.iter = 10,                  # One round of Harmony involves one clustering and one correction step
               #max.iter.cluster = 20,          # Maximum number of rounds to run clustering at each round of Harmony 
               early_stop = T,
               #dims.use = 30,
               )

#SeuratOBJ@reductions

message("Finishing Seurat-Harmony Integration - ", Sys.time())


## Run standard analysis: Calculate integrated.cca cell embeddings

# LogNormalize the count data present in the assay
SeuratOBJ <- NormalizeData(SeuratOBJ)

# Identifies features that are outliers on a 'mean variability plot'.
# vst method (default): First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess).
SeuratOBJ <- FindVariableFeatures(SeuratOBJ,
                                  selection.method = "vst") # First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression...

map(save_VFeatures, 'harmony_r')
# # Scale data
# # Run a PCA dimensionality reduction


SeuratOBJ <- FindNeighbors(SeuratOBJ, reduction = "harmony_r", dims = 1:30)
SeuratOBJ <- FindClusters(SeuratOBJ, resolution = 1)

SeuratOBJ <- RunUMAP(SeuratOBJ, dims = 1:30, reduction = "harmony_r")

SeuratOBJ@reductions$harmony_r

## Visualization
#reduct <- c('pca', 'umap.unintegrated', 'integrated.cca', 'harmony_r')
Map(plot_clust, reduct)



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

