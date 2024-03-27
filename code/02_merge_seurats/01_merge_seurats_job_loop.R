########################################################################
## Merge Seurat objects from an specific directory
##
## Authors. CSC
## Date. Feb 22, 2024
## Last.md: xxx
##
## Input: Seurat RDS Object generated with 01_preprocessing_GEX_ATAC.R
## Output:  New Seurat combined object: CCA and Harmony
##          Basic plots for reference after correction    
##
## NOTES: 30G free mem recommended for 20K cells
## For slurm env: runsrun --x11 --pty --partition=interactive bash
########################################################################

library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
#options(Seurat.object.assay.version = 'v5')    # To use new Seurat v5: Please run: options(Seurat.object.assay.version = 'v5')
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(here)

here::here()

# Check if processed_data directory exists, if not create it
if (!dir.exists(here("processed-data/02_merge_seurats/"))) {
    dir.create(here("processed-data/02_merge_seurats/"))
}
# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/02_merge_seurats/"))) {
    dir.create(here("plots/02_merge_seurats/"))
}
## Directory to save variable features 
dir <- file.path(here('processed-data/02_merge_seurats/csv_files/')) 
if (!dir.exists(dir)) dir.create(dir)





########################    Initials ########################  

## select the count-mtx to merge (raw or normalized data)
# count_mtx_type_label <- 'data_counts'      
# count_mtx_type_label <- 'norm_counts' 

## get args

args = commandArgs(trailingOnly=TRUE)
## read count mtx type (abs_counts and normalized_counts)
#count_mtx_type_label <- 'data_counts'      
count_mtx_type_label <- args[1]


## compose rds base name for merge data
if (count_mtx_type_label=='data_counts') { 
  s_sample <- 'seurat.combined.data_counts' 
} else { 
  s_sample <- 'seurat.combined.norm_counts' 
}

## read directory with seurat obj
path_directory_in <- here('processed-data/01_preprocessing_QC')
#filenames <- list.files(path_directory_in, pattern="*_Hb_KDM.rds")
all_rds <- paste0(path_directory_in, '/', list.files(path_directory_in, pattern="*_Hb_KDM.rds", recursive = TRUE))


## load pre-existing Seurat
get_seurat <- function(name) {

    #Ex. /dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/02_merge_seurats/seurat.combined.data_counts_PCA.rds
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
    png_name <- here('plots/02_merge_seurats', png_file)  
    ggsave(p1, filename = png_name, height = 5, width = 10)
    
    # visualize the two conditions side-by-side
    p1 <- DimPlot(sobj, 
                  reduction = reduct, split.by = "orig.ident")
    png_file <- paste0(f_name,'_dimplot_splitted.png')
    png_name <- here('plots/02_merge_seurats', png_file)  
    ggsave(p1, filename = png_name, height = 5, width = 10)
    
}


## function to save variable features before correction 
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






######## load pre-existing Seurat(s) objects not integrated / or Jump to load pre-existing Seurat integrated. ########

# SeuratOBJ <- get_seurat(here('processed-data/01_preprocessing_QC', 'S1_Hb_KDM.rds'))
# SeuratOBJ2 <- get_seurat(here('processed-data/01_preprocessing_QC', 'S2_Hb_KDM.rds'))


## Merge Seurat objects according with the `count type`
## NOTE: By default, merge() will combine the Seurat objects based on the raw count matrices, erasing any previously normalized and scaled data matrices. If you want to merge the normalized data matrices as well as the raw count matrices, simply pass merge.data = TRUE. This should be done if the same normalization approach was applied to all objects.

seurat_lst <- list()
seurat_name_lst <- list()

message('Starting to prepare seurats to merge ...')

for (rds_path in all_rds) {
  
  print(rds_path)
  SeuratOBJ <-readRDS(rds_path)
  
  if (count_mtx_type_label=='data_counts') {

    ## None additional step required

  } else {

    SeuratOBJ <- NormalizeData(SeuratOBJ)

  }
  
  if ( length(seurat_lst)>0 ) { seurat_lst <- append(seurat_lst, SeuratOBJ) } else { seurat_lst <- SeuratOBJ }
  if ( length(seurat_name_lst)>0 ) { seurat_name_lst <- append(seurat_name_lst, basename(rds_path)) } else { seurat_name_lst <- basename(rds_path) }
  
}

print(seurat_lst)
print(seurat_name_lst)

message('Seurats prepared')




## Merge object according with the given list(s)

# SeuratOBJ.combined <- merge(SeuratOBJ, y = SeuratOBJ2,
#                             add.cell.ids = c('S1_Hb', 'S2_Hb'),
#                             project = "Habenula", 
#                             merge.data = TRUE)     #  merge the normalized data matrices as well as the raw count matrices
# LayerData(SeuratOBJ.combined)[1:10, 1:15]
#pbmc.big <- merge(pbmc3k, y = c(pbmc4k, pbmc8k), add.cell.ids = c("3K", "4K", "8K"), project = "PBMC15K")


message('Merge completed!')









## verification of the integration

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

## Save integrated Seurat RDS Object

rds_name <- here('processed-data/02_merge_seurats', paste0(s_sample, '.rds'))
# .../seurat.combined.data_counts.rds"
saveRDS(SeuratOBJ.combined, file = rds_name)
message('Seurat combined saved in ', rds_name)   




############ Plots ############


## Violin plot with UMIs, Genes, ^MT and RIBO levels

p1 <- plot_violinQC(SeuratOBJ.combined, c("nCount_RNA", "nFeature_RNA", "percent.mt"), s_sample) 
#p1 <- VlnPlot(SeuratOBJ.combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "orig.ident") 
png_file <- paste0(s_sample, '_Vplots_GEX.png')
png_name <- here('plots/02_merge_seurats', png_file)  
ggsave(p1, filename = png_name, height = 5, width = 7)


## Plot Genes and UMIs by density per cell 

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
png_name <- here('plots/02_merge_seurats', png_file)  
ggsave(p1, filename = png_name, height = 4, width = 5)
message('UMI/Counts by MT plot saved!')  


## Gene distribution by cell

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
png_name <- here('plots/02_merge_seurats', png_file)  
ggsave(p1, filename = png_name, height = 4, width = 5)
message('UMI/Counts by MT plot saved!')  


## Correlation btw Genes/UMIs 

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
png_name <- here('plots/02_merge_seurats', png_file)  
ggsave(p1, filename = png_name, height = 4, width = 5)
message('UMI/Counts by MT plot saved!')  



############ Process PCAs

# NOTE timoast comment: https://github.com/satijalab/seurat/issues/3505
# I'd suggest doing QC and filtering cells on each object before running the integration.
# Running NormalizeData on the integrated assay will overwrite the integration results.

## Start from here IF yoi have load pre-existing combined Seurat(s)

if (count_mtx_type_label=='data_counts') { s_sample <- 'seurat.combined.data_counts' } else { s_sample <- 'seurat.combined.norm_counts' }
rds_name <- here('processed-data/02_merge_seurats', paste0(s_sample, '.rds'))
SeuratOBJ.combined <- get_seurat(rds_name)


##  verification of the integration

table(SeuratOBJ.combined$orig.ident)
# S1_Hb_KDM S2_Hb_KDM 
# 8178      9816 
head(colnames(SeuratOBJ.combined))
tail(colnames(SeuratOBJ.combined))
colnames(SeuratOBJ.combined@meta.data)


######### Perform analysis without integration

## split the RNA measurements into two layers one for each sample

SeuratOBJ.combined[["RNA"]] <- split(SeuratOBJ.combined[["RNA"]], f = SeuratOBJ.combined$orig.ident)
# Warning: Assay RNA changing from Assay to Assay5


## Run standard analysis: Calculate PCA cell embeddings

##  Perform logNormalize forthe count data present in the assay

SeuratOBJ.combined <- NormalizeData(SeuratOBJ.combined)


## Identifies features that are outliers on a 'mean variability plot'.
## vst method (default): First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess).

SeuratOBJ.combined <- FindVariableFeatures(SeuratOBJ.combined, selection.method = "vst") 
# Fits a line to the relationship of log(variance) and log(mean) using local polynomial regression...

## Save top 10,20,50 and 100 VF 
save_VFeatures(SeuratOBJ.combined, 'pca')


## Scale data. Perform “LogNormalize” method to the GEX for each cell by the total expression multiply by a scale factor (10,000 by default), and log-transforms the result

all.genes <- rownames(SeuratOBJ.combined)
SeuratOBJ.combined <- ScaleData(SeuratOBJ.combined, features = all.genes)


## Run a PCA dimensionality reduction

SeuratOBJ.combined <- RunPCA(SeuratOBJ.combined)
#SeuratOBJ.combined@reductions

p1 <- ElbowPlot(SeuratOBJ.combined)
png_file <- paste0(s_sample, '_PCAelbow.png')
png_name <- here('plots/02_merge_seurats', png_file)  
ggsave(p1, filename = png_name, height = 4, width = 5)


## Before data correction perform some visualizations for reference

## Compute nearest neighbor graph + SNN

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
png_name <- here('plots/02_merge_seurats', png_file)
ggsave(p1, filename = png_name) #, height = 5, width = 5

# Save RDS Object
rds_name <- here('processed-data/02_merge_seurats', paste0(s_sample, '_PCA.rds'))
# .../seurat.combined.data_counts_PCA.rds
saveRDS(SeuratOBJ.combined, file = rds_name)
message('Seurat combined saved in ', rds_name)   





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

