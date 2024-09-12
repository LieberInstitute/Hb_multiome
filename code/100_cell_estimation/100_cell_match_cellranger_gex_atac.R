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

# library(Seurat)                                 # 4.9.9.9045 2023-05-17 [1] Github (satijalab/seurat@7d1094c)
# #options(Seurat.object.assay.version = 'v5')    # To use new Seurat v5: Please run: options(Seurat.object.assay.version = 'v5')
# library(Signac)                                 # 1.9.0.9000 2023-05-08 [1] Github (stuart-lab/signac@cf31022)
 options(tidyverse.quiet = TRUE)
library(tidyverse)
library(readr)
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

### Barcode translation in Cell Rangers ###  

##    ATAC barcode: <path_to_cellranger-arc>/lib/python/atac/barcodes
##    GEX barcode: <path_to_cellranger-arc>/lib/python/cellranger/barcodes
##        The two sets of barcodes are associated by line number

## GEX barcodes
barcodes_GEX <- read.csv(here("/jhpce/shared/libd/core/cellranger_arc/2.0.2/cellranger-arc-2.0.2/lib/python/cellranger/barcodes", 
                                "737K-arc-v1.txt.gz"), header = FALSE)
# add row number as column
v_barcodes_GEX <- unlist(barcodes_GEX)
df_barcodes_GEX <- unlist(barcodes_GEX[1])
# barcodes_GEX$line <- 1:nrow(barcodes_GEX)
# tail(barcodes_GEX)
# colnames(barcodes_GEX)

# Split name column into firstname and last name


## ATAC barcodes
barcodes_ATAC <- read.csv(here("/jhpce/shared/libd/core/cellranger_arc/2.0.2/cellranger-arc-2.0.2/lib/python/atac/barcodes",
                               "737K-arc-v1.txt.gz"), header = FALSE)
#head(barcodes_ATAC)             


########################   Start   ########################  

## scan the sample's names

sample_args <- commandArgs(trailingOnly = TRUE)
# sample_name <- sample_args[1]
# For testing: 
#sample_args <- "4C_Hb_KDM,4A_Hb_KDM"
#sample_args <- "6C_Hb_KDM,6A_Hb_KDM"

sample_data <- unlist(strsplit(sample_args, ","))
sample_name_RNA <- sample_data[1]
sample_name_ATAC <- sample_data[2]

#sample_name <- sample_args[1]
message('Processing sample: ', sample_name_RNA)
b_get_filtered_GEX <- FALSE
b_get_GEX_plots <- FALSE


### Read barcodes RNA - only detected cell-associated barcodes

# feature.RNA.loc_Dir <- here(cellrangerDir_GEX, sample_name_RNA, "outs", "filtered_feature_bc_matrix", "features.tsv.gz")
# df_sample_features_RNA <- read.csv(feature.RNA.loc_Dir, header = FALSE, sep = "\t")
# head(df_sample_features_RNA)

barcode.RNA.loc_Dir <- here(cellrangerDir_GEX, sample_name_RNA, "outs", 
                            "filtered_feature_bc_matrix", "barcodes.tsv.gz")
df_sample_barcodes_RNA <- read.csv(barcode.RNA.loc_Dir, header = FALSE)
head(df_sample_barcodes_RNA)
# remove GEM-Bed number `-1` and creates a vector
v_sample_barcodes_RNA <- str_remove(unlist(df_sample_barcodes_RNA), "-1")
v_sample_barcodes_RNA <- unlist(as.data.frame(v_sample_barcodes_RNA))

# # Validate and test: 
# head(v_barcodes_GEX)
# length(v_barcodes_GEX)
# head(v_sample_barcodes_RNA)
# length(v_sample_barcodes_RNA)
# v_sample_barcodes_RNA[5] <- v_barcodes_GEX[10]  #AAACAGCCAAACTGCC
# v_sample_barcodes_RNA[10] <- v_barcodes_GEX[100] #AAACAGCCAATTATGC

matching_cells_gex <- pmatch(v_sample_barcodes_RNA, v_barcodes_GEX)
message("Translated barcodes from cellranger-arc GEX: ", length(matching_cells_gex))

# head(matching_cells_gex, n=10)
# v_sample_barcodes_RNA[10]==v_barcodes_GEX[100]




### Read barcodes ATAC - - only detected cell-associated barcodes

barcode.ATAC.loc_Dir <- here(cellrangerDir_ATAC, sample_name_ATAC, "outs", 
                             "filtered_peak_bc_matrix", "barcodes.tsv")
df_sample_barcodes_ATAC <- readr::read_tsv(barcode.ATAC.loc_Dir, 
                                           trim_ws = TRUE , 
                                           col_names = FALSE, show_col_types = FALSE)
head(df_sample_barcodes_ATAC)
v_sample_barcodes_atac <- str_remove(unlist(df_sample_barcodes_ATAC), "-1")
#v_sample_barcodes_atac[1:10]

## For testing equalities 
#v_sample_barcodes_atac[2] <- v_sample_barcodes_RNA[1]
#identical(v_sample_barcodes_atac[1], v_sample_barcodes_RNA[1]) # False
#identical(v_sample_barcodes_atac[2], v_sample_barcodes_RNA[1]) # True

message("Barcode lenghts (RNA/ATAC)")
length(v_sample_barcodes_RNA)
length(v_sample_barcodes_atac)

match_cells <- intersect(v_sample_barcodes_RNA, v_sample_barcodes_atac)
number_match_cells <- length(match_cells)
message("Matching cells: ", number_match_cells)

## Cells matched
write.csv(match_cells, row.names = FALSE, quote = FALSE, 
          here(processedDir, paste0(sample_name_RNA, "_", sample_name_ATAC, "_matching_cells")))

## number and percentage
percent_rna <- (number_match_cells * 100) / length(v_sample_barcodes_RNA)
percent_atac <- (number_match_cells * 100) / length(v_sample_barcodes_atac)
cvs_percents <- "Sample,Total_number_rna,Total_number_atac,Total_match,Percent_rna,Percent_atac\n"
cvs_percents <- paste0(cvs_percents,  paste0(paste0(sample_name_RNA, "-", sample_name_ATAC), ",",
                      length(v_sample_barcodes_RNA), ",", length(v_sample_barcodes_atac), ",",
                      number_match_cells, ",", percent_rna, percent_atac, "\n"))
cat(cvs_percents)
write.csv(cvs_percents, row.names = FALSE, quote = FALSE, 
          here(processedDir, paste0(sample_name_RNA, "_", sample_name_ATAC, "_percent_matching_cells")))


message(" Processed ", sample_name_RNA, " and ", sample_name_ATAC)


# ### Reads barcodes from seurat object

# filtered_barcode_path <- here(cellrangerDir_GEX, sample_name_RNA, "outs", "filtered_feature_bc_matrix.h5")

# rna_counts <- Read10X_h5(filtered_barcode_path)             
# head(rna_counts, n = 3)
# SeuratOBJ = CreateSeuratObject(counts = rna_counts)
# message('Seurat object created successfully!')
# head(SeuratOBJ)
# # removes the "-1" if all cell names contain it
# SeuratOBJ <- RenameCells(SeuratOBJ, new.names = str_remove(Cells(x = SeuratOsBJ), "-1"))
# df_cells <- as.data.frame(SeuratOBJ@meta.data, row.names = NULL)
# v_sample_barcodes2 <- rownames(df_cells)
# #sample_barcodes2 <- data.table(sample_barcodes2)
# identical(v_sample_barcodes, v_sample_barcodes2)


############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()

