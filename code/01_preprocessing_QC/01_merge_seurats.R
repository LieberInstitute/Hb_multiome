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

#source(here("code/functions_custom", "remote_file_caller.R"))       # Call functions to read paths
#source(here("code/functions_custom", "remote_seurat_functions.R"))  # Call functions to create and handle Seurat object
#source(here("code/functions_custom", "remote_signac_functions.R"))  # Call functions to create Signac object
source(here("code/functions_custom", "remote_plot_functions.R"))    # Call to plot GEX assay
source(here("code/functions_custom", "remote_filtering_functions.R"))   # Call functions to subset the Seurat object


########################    Initials ########################  

# seurats_lst <- c('S1_Hb_KDM', 'S2_Hb_KDM')
# 
# bfirst <- FALSE
# for (x in seurats_lst) {
#     if (!bfirst) {
#         SeuratOBJ.combined <- readRDS(here('processed-data/01_preprocessing_QC', paste0(x, '.rds')))
#         bfirst <- TRUE
#     } else {
#         SeuratOBJ <- readRDS(here('processed-data/01_preprocessing_QC', paste0(x, '.rds')))
#         #SeuratOBJ@meta.data
#         #str(SeuratOBJ)
#         print(SeuratOBJ)
#         
#         SeuratOBJ.combined <- merge(SeuratOBJ.combined, y = SeuratOBJ,
#                                     add.cell.ids = c(s_sample, 'PBMC'),
#                                     project = "PBMC-FFB")
#     }
#     
#     table(SeuratOBJ.combined$orig.ident)
#     
# }


# load pre-existing seurat objects
SeuratOBJ <- readRDS(here('processed-data/01_preprocessing_QC', 'S1_Hb_KDM.rds'))

SeuratOBJ2 <- readRDS(here('processed-data/01_preprocessing_QC', 'S2_Hb_KDM.rds'))

table(SeuratOBJ$orig.ident)

# select the count-mtx to merge (raw or normalized data)
count_mtx_type <- 'raw_counts'      
#count_mtx_type <- 'normalized'      

# Merge Seurat objects according with the `count_mtx_type`
if (count_mtx_type=='raw_counts') {
    
    SeuratOBJ.combined <- merge(SeuratOBJ, y = SeuratOBJ2,
                                add.cell.ids = c('S1_Hb_KDM', 'S2_Hb_KDM'),
                                project = "Habenula")
    LayerData(SeuratOBJ.combined)[1:10, 1:15]
    
} else {
    
    SeuratOBJ <- NormalizeData(SeuratOBJ)
    SeuratOBJ2 <- NormalizeData(SeuratOBJ2)
    SeuratOBJ.combined <- merge(SeuratOBJ, y = SeuratOBJ2,
                                  add.cell.ids = c('S1_Hb_KDM', 'S2_Hb_KDM'),
                                  project = "Habenula", 
                                  merge.data = TRUE)     #  merge the normalized data matrices as well as the raw count matrices
    LayerData(SeuratOBJ.combined)[1:10, 1:15]
}


#
#pbmc.big <- merge(pbmc3k, y = c(pbmc4k, pbmc8k), add.cell.ids = c("3K", "4K", "8K"), project = "PBMC15K")
message('Merge complete!')

# verification of the integration
table(SeuratOBJ.combined$orig.ident)
head(colnames(SeuratOBJ.combined))
tail(colnames(SeuratOBJ.combined))
unique(sapply(X = strsplit(colnames(SeuratOBJ.combined), split = "_"), FUN = "[", 1))

head(SeuratOBJ.combined)

# basic stats gor GEX
p1_ATAC <- VlnPlot(SeuratOBJ.combined, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "orig.ident") 
p1_ATAC

# Save RDS Object
rds_name <- here('processed-data/01_preprocessing_QC', 'seurat.combined_GEX.rds')
saveRDS(SeuratOBJ.combined, file = rds_name)
message('Seurat combined saved in ', rds_name)   


############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()

