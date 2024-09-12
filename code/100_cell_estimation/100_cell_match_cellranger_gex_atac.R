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
v_barcodes_GEX <- unlist(barcodes_GEX)
# barcodes_GEX$line <- 1:nrow(barcodes_GEX)

## ATAC barcodes
barcodes_ATAC <- read.csv(here("/jhpce/shared/libd/core/cellranger_arc/2.0.2/cellranger-arc-2.0.2/lib/python/atac/barcodes",
                               "737K-arc-v1.txt.gz"), header = FALSE)
v_barcodes_ATAC <- unlist(barcodes_ATAC)


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

message('Processing sample: ', sample_name_RNA)


### Read barcodes RNA - only detected cell-associated barcodes

# feature.RNA.loc_Dir <- here(cellrangerDir_GEX, sample_name_RNA, "outs", "filtered_feature_bc_matrix", "features.tsv.gz")
# df_sample_features_RNA <- read.csv(feature.RNA.loc_Dir, header = FALSE, sep = "\t")
# head(df_sample_features_RNA)

barcode.RNA.loc_Dir <- here(cellrangerDir_GEX, sample_name_RNA, "outs", 
                            "filtered_feature_bc_matrix", "barcodes.tsv.gz")
v_sample_bc_gex <- read.csv(barcode.RNA.loc_Dir, header = FALSE)
head(v_sample_bc_gex)
# remove GEM-Bed number `-1` and creates a vector
v_sample_bc_gex <- str_remove(unlist(v_sample_bc_gex), "-1")
df_gex <- as.data.frame(v_sample_bc_gex) # to add line-number corresponding barcode from cellranger-arc
v_sample_bc_gex <- unlist(as.data.frame(v_sample_bc_gex))

# # Validate and test: 
# head(v_barcodes_GEX)
# length(v_barcodes_GEX)
# head(v_sample_bc_gex)
# length(v_sample_bc_gex)
# v_sample_bc_gex[5] <- v_barcodes_GEX[10]  #AAACAGCCAAACTGCC
# v_sample_bc_gex[10] <- v_barcodes_GEX[100] #AAACAGCCAATTATGC

translated_bc_arc_gex <- pmatch(v_sample_bc_gex, v_barcodes_GEX)
message("Translated barcodes from cellranger-arc GEX: ", length(translated_bc_arc_gex))

# head(translated_bc_arc_gex, n=10)
# v_sample_bc_gex[10]==v_barcodes_GEX[100]

## Add line-number corresponding barcode from cellranger-arc
df_gex$ID <- translated_bc_arc_gex
head(df_gex)


### Read barcodes ATAC - - only detected cell-associated barcodes

barcode.ATAC.loc_Dir <- here(cellrangerDir_ATAC, sample_name_ATAC, "outs", 
                             "filtered_peak_bc_matrix", "barcodes.tsv")
df_sample_bc_atac <- readr::read_tsv(barcode.ATAC.loc_Dir, 
                                           trim_ws = TRUE , 
                                           col_names = FALSE, show_col_types = FALSE)
head(df_sample_bc_atac)
v_sample_bc_atac <- str_remove(unlist(df_sample_bc_atac), "-1")
df_atac <- as.data.frame(v_sample_bc_atac) # to add line-number corresponding barcode from cellranger-arc
v_sample_bc_atac <- unlist(as.data.frame(v_sample_bc_atac))

translated_bc_arc_atac <- pmatch(v_sample_bc_atac, v_barcodes_ATAC)
message("Translated barcodes from cellranger-arc ATAC: ", length(translated_bc_arc_atac))
#head(translated_bc_arc_atac)
# Testing
v_sample_bc_atac[1] == v_barcodes_ATAC[473848]
v_sample_bc_atac[5] == v_barcodes_ATAC[550741]

## match translated number barcodes between GEX and ATAC assays alone
match_cells <- intersect(translated_bc_arc_gex, translated_bc_arc_atac)
number_match_cells <- length(match_cells)
message("Matching cells: ", number_match_cells)

## Add line-number corresponding barcode from cellranger-arc
df_atac$ID <- translated_bc_arc_atac


## Join table with translated barcodes matching btw both assays
df_matching_barcodes <- inner_join(df_gex, df_atac, by = "ID")

# ## Additional validation
## For ATAC
sum_atac <- sum(translated_bc_arc_atac %in% df_matching_barcodes$ID)
sum_atac_valid_bc <- sum(sum_atac, na.rm = TRUE)
sum_atac_NO_valid_bc <- length(translated_bc_arc_atac) - sum_atac_valid_bc

sum_gex <- sum(translated_bc_arc_gex %in% df_matching_barcodes$ID)
sum_gex_valid_bc <- sum(sum_atac, na.rm = TRUE)
sum_gex_NO_valid_bc <- length(translated_bc_arc_gex) - sum_gex_valid_bc

if (sum_gex_valid_bc == sum_atac_valid_bc) {
  ## Cells matched
  write.csv(df_matching_barcodes, row.names = TRUE, quote = FALSE, 
            here(processedDir, paste0(sample_name_RNA, "_", sample_name_ATAC, "_matching_cells.csv")))
  message("Matching ARC translated barcodes saved!")
} else {
  message("Error translating barcodes!")
}
  
(nrow(df_matching_barcodes) == number_match_cells & length(rowSums(is.na(df_matching_barcodes)) == 0))

## Calculate percentages
if ((sum_atac_valid_bc + sum_atac_NO_valid_bc) == length(translated_bc_arc_atac)) {
  percent_atac <- (sum_atac_valid_bc * 100) / length(v_sample_bc_atac)  
  message("Percentage of matching ATAC translated barcodes: ", percent_atac)
} else {
  message("Translated barcodes from ATAC does not match")
}
 
if ((sum_gex_valid_bc + sum_gex_NO_valid_bc) == length(translated_bc_arc_gex)) {
  percent_gex <- (sum_gex_valid_bc * 100) / length(v_sample_bc_gex)  
  message("Percentage of matching GEX translated barcodes: ", percent_gex)
} else {
  message("Translated barcodes from GEX does not match")
}

## Saving basic stats 
tsv_header <- "Sample_id\tsample_gex_bc\tvalid_bc_gex_arc\tp_valid_bc_gex_arc\tnot_valid_bc_gex_arc\tp_not_valid_bc_gex_arc\tvalid_bc_atac_arc\tp_valid_bc_atac_arc\tnot_valid_bc_atac_arc\tp_not_valid_bc_atac_arc\n"
s_name <- paste0(paste0(sample_name_RNA, "-", sample_name_ATAC), "\t")
gex_line <- paste0(length(translated_bc_arc_gex), "\t", sum_gex_valid_bc, "\t", percent_gex, "\t", sum_gex_NO_valid_bc, "\t", (100-percent_gex), "\t")  
atac_line <- paste0(length(translated_bc_arc_atac), "\t", sum_atac_valid_bc, "\t", percent_atac, "\t", sum_atac_NO_valid_bc, "\t", (100-percent_atac))  
tsv_body <- paste0(tsv_header, s_name, gex_line, atac_line)
cat(tsv_body)
write.table(tsv_body, here(processedDir, paste0(sample_name_RNA, "_", sample_name_ATAC, "_matching_cells_stats.tsv")),
            quote=FALSE, sep='\t', col.names = NA)

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

