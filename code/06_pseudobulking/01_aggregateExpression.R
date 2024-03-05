########################################################################
## Aggregateexpression for pseudobulking
## Authors. CSC
## Date. March 5th, 2024
## Last.Adaptation: xxx
##
## Input: Seurat integrated object with samples S1 and S2 after CCA or Harmony correction
## Output:  
##
## NOTES:
## For slurm env: runsrun --x11 --pty --partition=interactive bash
########################################################################

library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
#options(Seurat.object.assay.version = 'v5')    # To use new Seurat v5: Please run: options(Seurat.object.assay.version = 'v5')
#library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
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
if (!dir.exists(here("processed-data/06_pseudobulking/"))) {
    dir.create(here("processed-data/06_pseudobulking/"))
}
# Check if plot directory exists, if not create it
if (!dir.exists(here("plots/06_pseudobulking/"))) {
    dir.create(here("plots/06_pseudobulking/"))
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



# load pre-existing seurat objects
if (count_mtx_type=='data_counts') { s_sample <- 'seurat.combined.data_counts_PCA_CCA_Harmony.rds' } else { s_sample <- 'seurat.combined.norm_counts_PCA_CCA_Harmony.rds' }
rds_name <- here('processed-data/01_preprocessing_QC', paste0(s_sample, '.rds'))
SeuratOBJ <- get_seurat(rds_name)

# verification of the integration
table(SeuratOBJ$orig.ident)
# S1_Hb_KDM S2_Hb_KDM 
# 8178      9816 
head(colnames(SeuratOBJ))
tail(colnames(SeuratOBJ))
colnames(SeuratOBJ@meta.data)


######### Perform aggregation 

## https://satijalab.org/seurat/reference/aggregateexpression for pseudobulking
##     and then use the DoHeatmap() function they have https://satijalab.org/seurat/reference/doheatmap
##.    https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/pseudobulk-expression.html

AggregateExpression(
  object,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  margin = 1,
  verbose = TRUE,
)


######### Plot Heatmap with DoHeatmep()

# https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/pseudobulk-expression.html

# https://alexyermanos.github.io/Platypus/reference/GEX_cluster_genes_heatmap.html#:~:text=Produces%20a%20heatmap%20displaying%20the,edited%20using%20typical%20ggplot%20interactions.


DoHeatmap(
  object,
  features = NULL,
  cells = NULL,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  vjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)


png_file <- paste0(s_sample, '_integrated.cca_pca_heatmap.png')
png_name <- here('plots/04_merge_seurats', png_file)
ggsave(p1, filename = png_name, height = 5, width = 10)


rds_name <- here('processed-data/04_merge_seurats', paste0(s_sample, '_PCA_CCA_Harmony.rds'))
# .../seurat.combined.data_counts_PCA_CCA_Harmony.rds
saveRDS(SeuratOBJ, file = rds_name)
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

