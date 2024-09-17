########################################################################
## Create csv file with matching ARC barcodes for GEX and ATAC analyzed by separate.
## Create tsv file with basic stats derived from this analysis.
## Directions: https://kb.10xgenomics.com/hc/en-us/articles/360049105612-Barcode-translation-in-Cell-Ranger-ARC

## Authors. CSC
## Date. Sep 10th, 2024
## NOTES: 
## GEX was processed with cellranger-count from the multiome datasets (GEX only)
## ATAC was processed with cellranger-atac from the multiome datasets (ATAC only)
##
## For slurm env: $srun --pty --mem=15GB --x11 bash
#############################################################

library(tidyverse)
library(readr)
library(here)


########################    Initials ########################  

### working directories
here::here()
main_dir_name <- "100_cell_match_cellranger_gex_atac"
cellrangerDir_GEX <- here("processed-data", "cellrangerGEX")
cellrangerDir_ATAC <- here("processed-data", "cellrangerATAC")
cellrangerDir_reanalize <- here("processed-data", main_dir_name, "reanalize_files")
functionsDir <- here("code", "functions_custom")
processedDir <- here("processed-data", main_dir_name)

# Check Directories
if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(cellrangerDir_reanalize)) { dir.create(cellrangerDir_reanalize) }

### Barcode translation references in Cell Ranger ARC lives here
##    ATAC barcode: <path_to_cellranger-arc>/lib/python/atac/barcodes
##    GEX barcode: <path_to_cellranger-arc>/lib/python/cellranger/barcodes
##        The two sets of barcodes are associated by line number

## Load GEX and ATAC barcodes
barcodes_GEX <- read.csv(here("/jhpce/shared/libd/core/cellranger_arc/2.0.2/cellranger-arc-2.0.2/lib/python/cellranger/barcodes", 
                                "737K-arc-v1.txt.gz"), header = FALSE)
v_barcodes_GEX <- unlist(barcodes_GEX)

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
#sample_args <- "3C_Hb_KDM,3A_Hb_KDM"

sample_data <- unlist(strsplit(sample_args, ","))
sample_name_RNA <- sample_data[1]
sample_name_ATAC <- sample_data[2]

message('Processing sample: ', sample_name_RNA)


### Translate RNA barcodes - only detected cell-associated barcodes

barcode.RNA.loc_Dir <- here(cellrangerDir_GEX, sample_name_RNA, "outs", 
                            "filtered_feature_bc_matrix", "barcodes.tsv.gz")
v_sample_bc_gex <- read.csv(barcode.RNA.loc_Dir, header = FALSE)
#head(v_sample_bc_gex)
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

## retrieve line-number of matching vector in the GEX cellranger-arc reference
translated_bc_arc_gex <- pmatch(v_sample_bc_gex, v_barcodes_GEX)
message("Translated barcodes from cellranger-arc GEX: ", length(translated_bc_arc_gex))

## Add line-number corresponding translated ARC barcode
df_gex$ID <- translated_bc_arc_gex

### Translate ATAC barcodes - only detected cell-associated barcodes

barcode.ATAC.loc_Dir <- here(cellrangerDir_ATAC, sample_name_ATAC, "outs", 
                             "filtered_peak_bc_matrix", "barcodes.tsv")
df_sample_bc_atac <- readr::read_tsv(barcode.ATAC.loc_Dir, trim_ws = TRUE, 
                                           col_names = FALSE, show_col_types = FALSE)
#head(df_sample_bc_atac)
v_sample_bc_atac <- str_remove(unlist(df_sample_bc_atac), "-1")
df_atac <- as.data.frame(v_sample_bc_atac) # to add line-number corresponding barcode from cellranger-arc
v_sample_bc_atac <- unlist(as.data.frame(v_sample_bc_atac))

translated_bc_arc_atac <- pmatch(v_sample_bc_atac, v_barcodes_ATAC)
message("Translated barcodes from cellranger-arc ATAC: ", length(translated_bc_arc_atac))

## match translated barcodes between GEX and ATAC assays alone
match_cells <- intersect(translated_bc_arc_gex, translated_bc_arc_atac)
number_match_cells <- length(match_cells)
message("Matching cells: ", number_match_cells)

# head(sort(translated_bc_arc_gex), n=10)
# head(sort(translated_bc_arc_atac), n=10)
# head(match_cells, n=10)

## Add line-number corresponding translated ARC barcode
df_atac$ID <- translated_bc_arc_atac

## Join table with translated barcodes matching both assays
df_matching_barcodes <- inner_join((df_gex |> arrange(ID)), (df_atac |> arrange(ID)), by = "ID")
# head(df_matching_barcodes)

## Compute validations and prepare some vectors to compute basic stats
## For ATAC
sum_atac <- sum(translated_bc_arc_atac %in% df_matching_barcodes$ID)
sum_atac_valid_bc <- sum(sum_atac, na.rm = TRUE)
sum_atac_NO_valid_bc <- length(translated_bc_arc_atac) - sum_atac_valid_bc
## For GEX
sum_gex <- sum(translated_bc_arc_gex %in% df_matching_barcodes$ID)
sum_gex_valid_bc <- sum(sum_gex, na.rm = TRUE)
sum_gex_NO_valid_bc <- length(translated_bc_arc_gex) - sum_gex_valid_bc

if (sum_gex_valid_bc == sum_atac_valid_bc) {
  ## Matching cells between cellranger-count and cellranger-atac 
  ## In csv file: col1: gex barcode, col2: arc line number, col3: atac barcode
  names(df_matching_barcodes) <- c("v_sample_bc_gex", "cellrangerARC_line_ID", "v_sample_bc_atac")
  write.csv(df_matching_barcodes, row.names = TRUE, quote = FALSE, 
            here(processedDir, paste0(sample_name_RNA, "_", sample_name_ATAC, "_matching_cells.csv")))
  message("Matching ARC translated barcodes saved!")

  ## Barcodes from `cellranger ARC` matching between cellranger-count and cellranger-atac 
  ## This selection of barcodes matching is for running with `cellranger-arc reanalyze` pipeline
  ARC_barcodes_matching <- v_barcodes_GEX[c(match_cells)]
  if (length(ARC_barcodes_matching) == sum_gex_valid_bc) {
    write.csv(as.data.frame(ARC_barcodes_matching), row.names = FALSE, quote = FALSE, 
              here(cellrangerDir_reanalize, paste0(sample_name_RNA, "_matching_ARC_barcodes.csv")))
    message("Matching ARC translated barcodes saved!")
  } else {
    message("Error saving cellranger-ARC barcodes!")
  }
} else {
  message("Error translating barcodes!")
}
  
## Calculate percentages
if (((sum_atac_valid_bc + sum_atac_NO_valid_bc) == length(translated_bc_arc_atac)) &  
  ((sum_gex_valid_bc + sum_gex_NO_valid_bc) == length(translated_bc_arc_gex))) {
    percent_atac <- (sum_atac_valid_bc * 100) / length(v_sample_bc_atac)  
    message("Percentage of matching ATAC translated barcodes: ", percent_atac)
    percent_gex <- (sum_gex_valid_bc * 100) / length(v_sample_bc_gex)  
    message("Percentage of matching GEX translated barcodes: ", percent_gex)
} else {
    message("Translated GEX and ATAC ARC barcodes do not match")
}

## Save the tsv file with basic stats 
tsv_header <- "Sample_id\t" 
tsv_header <- paste0(tsv_header, "translated_gex_bc\t valid_bc_gex_arc\t p_valid_bc_gex_arc\t not_valid_bc_gex_arc\t p_not_valid_bc_gex_arc\t")
tsv_header <- paste0(tsv_header, "translated_atac_bc\t valid_bc_atac_arc\t p_valid_bc_atac_arc\t not_valid_bc_atac_arc\t p_not_valid_bc_atac_arc\n")
s_name <- paste0(paste0(sample_name_RNA, "-", sample_name_ATAC), "\t")
gex_line <- paste0(length(translated_bc_arc_gex), "\t", sum_gex_valid_bc, "\t", percent_gex, "\t", sum_gex_NO_valid_bc, "\t", (100-percent_gex), "\t")  
atac_line <- paste0(length(translated_bc_arc_atac), "\t", sum_atac_valid_bc, "\t", percent_atac, "\t", sum_atac_NO_valid_bc, "\t", (100-percent_atac))  
tsv_body <- paste0(tsv_header, s_name, gex_line, atac_line)
cat(tsv_body)
write.table(tsv_body, here(processedDir, paste0(sample_name_RNA, "_", sample_name_ATAC, "_matching_cells_stats.tsv")),
            quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)

message(" Processed ", sample_name_RNA, " and ", sample_name_ATAC)




############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()

