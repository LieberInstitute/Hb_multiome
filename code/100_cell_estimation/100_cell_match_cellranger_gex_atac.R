########################################################################
## Gets the matching barcodes of multi-ome cell estimation between cellranger-count and cellranger-atac 
## Directions: https://kb.10xgenomics.com/hc/en-us/articles/360049105612-Barcode-translation-in-Cell-Ranger-ARC

## Authors. CSC
## Date. Sep 10th, 2024
## Last Update: xx
##
## NOTES: 
## GEX was processed with cellranger-count from the multiome datasets (GEX only)
## ATAC was processed with cellranger-atac from the multiome datasets (ATAC only)
##
## For slurm env: $srun --pty --mem=40GB --x11 bash
########################################################################

library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
#options(Seurat.object.assay.version = 'v5')    # To use new Seurat v5: Please run: options(Seurat.object.assay.version = 'v5')
library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(here)


########################    Initials ########################  

### working directories
here::here()
main_dir_name <- "100_cell_match_cellranger_gex_atac"
cellrangerDir_GEX <- here("processed-data", "cellrangerGEX")
cellrangerDir_ATAC <- here("processed-data", "cellrangerATAC")
functionsDir <- here("code", "functions_custom")
processedDir <- here("processed-data", main_dir_name)
#plotDir <- here("plots", main_dir_name)

# Check processed_data and plot directories exists
if (!dir.exists(processedDir)) { dir.create(processedDir) }
#if (!dir.exists(plotDir)) { dir.create(plotDir) }
source(here(functionsDir, "remote_seurat_functions.R"))  # Call functions to create and handle Seurat object

### Barcode translation in Cell Rangers ###  

##    ATAC barcode: <path_to_cellranger-arc>/lib/python/atac/barcodes
##    GEX barcode: <path_to_cellranger-arc>/lib/python/cellranger/barcodes
##        The two sets of barcodes are associated by line number

## GEX barcodes

barcodes_GEX <- read.csv(here("/jhpce/shared/libd/core/cellranger_arc/2.0.2/cellranger-arc-2.0.2/lib/python/cellranger/barcodes", 
                                "737K-arc-v1.txt.gz"), header = FALSE)
head(barcodes_GEX)

## ATAC barcodes

barcodes_ATAC <- read.csv(here("/jhpce/shared/libd/core/cellranger_arc/2.0.2/cellranger-arc-2.0.2/lib/python/atac/barcodes",
                               "737K-arc-v1.txt.gz"), header = FALSE)
head(barcodes_ATAC)             


########################   Start   ########################  

## scan the sample's names

sample_args <- commandArgs(trailingOnly = TRUE)
# sample_name <- sample_args[1]
# For testing: sample_args <- "4C_Hb_KDM,4A_Hb_KDM"

sample_data <- unlist(strsplit(sample_args, ","))
sample_name_RNA <- sample_data[1]
sample_name_ATAC <- sample_data[2]

#sample_name <- sample_args[1]
message('Processing sample: ', sample_name_RNA)
b_get_filtered_GEX <- FALSE
b_get_GEX_plots <- FALSE


### Read barcodes RNA 

barcode.RNA.loc_Dir <- here(cellrangerDir_GEX, sample_name_RNA, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")

### Pull barcodes from barcode matrix
df_sample_barcodes_RNA <- read.csv(barcode.RNA.loc_Dir, header = FALSE)
head(df_sample_barcodes_RNA)
v_sample_barcodes <- str_remove(unlist(df_sample_barcodes_RNA), "-1")
v_sample_barcodes[1:10]
# df_sample_barcodes_RNA <- data.frame(v_sample_barcodes)
# head(df_sample_barcodes_RNA, n=3)


# ### Reads barcodes from seurat object

# filtered_barcode_path <- here(cellrangerDir_GEX, sample_name_RNA, "outs", "filtered_feature_bc_matrix.h5")

# rna_counts <- Read10X_h5(filtered_barcode_path)             
# head(rna_counts, n = 3)
# SeuratOBJ = CreateSeuratObject(counts = rna_counts)
# message('Seurat object created successfully!')
# head(SeuratOBJ)
# # removes the "-1" if all cell names contain it
# SeuratOBJ <- RenameCells(SeuratOBJ, new.names = str_remove(Cells(x = SeuratOBJ), "-1"))
# df_cells <- as.data.frame(SeuratOBJ@meta.data, row.names = NULL)
# v_sample_barcodes2 <- rownames(df_cells)
# #sample_barcodes2 <- data.table(sample_barcodes2)
# identical(v_sample_barcodes, v_sample_barcodes2)

head(v_sample_barcodes)
head(unlist(barcodes_GEX))


### Read barcodes ATAC 

barcode.ATAC.loc_Dir <- here(cellrangerDir_ATAC, sample_name_ATAC, "outs", "filtered_peak_bc_matrix", "barcodes.tsv")
df_sample_barcodes_ATAC <- read.csv(barcode.ATAC.loc_Dir, header = FALSE)
head(df_sample_barcodes_ATAC)
v_sample_barcodes_atac <- str_remove(unlist(df_sample_barcodes_ATAC), "-1")
v_sample_barcodes_atac[1:10]
# df_sample_barcodes <- data.frame(v_sample_barcodes)
# head(df_sample_barcodes, n=3)











### Add additional meta-data: Chr-Mitochondrial levels

SeuratOBJ$log10GenesPerUMI <- log10(SeuratOBJ$nFeature_RNA) / log10(SeuratOBJ$nCount_RNA)
SeuratOBJ[["percent.mt"]] <- PercentageFeatureSet(SeuratOBJ, pattern = "^MT-")
SeuratOBJ[["percent.ribo"]] <- PercentageFeatureSet(SeuratOBJ, pattern = "^RP[LS]")
SeuratOBJ[["MTRatio"]] <- SeuratOBJ$percent.mt / 100 
message('Mitochondrial and Ribosomal percentage levels added')


########  ################################################# ######## 
########        2. Get Visualizations for the GEX           ######## 
########  ################################################# ######## 

# Build UMI, Genes, MITO and RIBO violin plots, number of cells per sample, UMI/transcripts per cell plots, and 
#   Distribution of genes per cell histogram

## Before QCed data, plot GEX basic quality controls
if (b_get_GEX_plots) {
    base_name <- paste0(sample_name, '_None_QC')
    plot_GEX_QCs(SeuratOBJ, base_name, TRUE)
    ## Calculate basic interquartile range for basic GEX stats
    # table_descriptive_stats_GEX(SeuratOBJ, base_name, s_tissue)
}


########  ################################################# ######## 
##        3. Preprocess Seurat Object based on the Mitochondrial percentage in the GEX assay
##              METHOD 1: (M1p)  calculate cut-off based on probabilities
##              METHOD 2: (M2sd) calculate cut-off at +/-2 Standard Deviations (95%)
########  ################################################# ######## 


# Create a second Seurat filtered by custom method of 1) probabilities or 2) SD
if (b_get_filtered_GEX) {
    # Preferred method 1. filter by probabilities
    
    # Validate or Load Seurat object from disk if available
    # Note. GEX assay need to be active

    if (i_filtering_method==1) {
        # filter by probabilities: M1
        SeuratOBJ.filtered <- get_seurat_GEX_filteringM1p(SeuratOBJ, FALSE)
    } else {    
        # filter by SD: M2
        SeuratOBJ.filtered <- get_seurat_GEX_filteringM2sd(SeuratOBJ, FALSE)
    }
    # Rename origin ident to identify the filtered Seurat object
    base_name <- paste0(levels(SeuratOBJ$`orig.ident`[1]), '_M', i_filtering_method)
    SeuratOBJ.filtered$orig.ident <- base_name
    #SeuratOBJ.filtered$`orig.ident`[1]

    plot_GEX_QCs(SeuratOBJ.filtered, base_name, TRUE)
    # Calculate and save some descriptive stats for further analysis
    table_descriptive_stats_GEX(SeuratOBJ.filtered, base_name, s_tissue)
    
    # Create a list of Seurat Objects to attach ATAC assay and process by QC
    lst_seurats <- list(SeuratOBJ, SeuratOBJ.filtered)
    
    message('Seurat filtered by MT level successfully!')

    rds_name <- here(processedDir, paste0(sample_name,'filteredM',i_filtering_method,'.rds'))
    saveRDS(SeuratOBJ.filtered, file = rds_name)
    message('Saving Seurat filtered.')
    
} else {
    
    # Only one Seurat object available
    lst_seurats <- list(SeuratOBJ)    
    rds_name <- here(processedDir, paste0(sample_name,'.rds'))
    saveRDS(SeuratOBJ, file = rds_name)
    message('Saving Seurat none filtered.')
    
}




############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()

