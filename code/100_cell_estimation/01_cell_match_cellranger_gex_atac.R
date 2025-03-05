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
library(Seurat)
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
barcodes_GEX_ARC <- read.csv(here("/jhpce/shared/libd/core/cellranger_arc/2.0.2/cellranger-arc-2.0.2/lib/python/cellranger/barcodes", 
                                "737K-arc-v1.txt.gz"), header = FALSE)
v_barcodes_GEX_ARC <- unlist(barcodes_GEX_ARC)

barcodes_ATAC_ARC <- read.csv(here("/jhpce/shared/libd/core/cellranger_arc/2.0.2/cellranger-arc-2.0.2/lib/python/atac/barcodes",
                               "737K-arc-v1.txt.gz"), header = FALSE)
v_barcodes_ATAC_ARC <- unlist(barcodes_ATAC_ARC)


########################   Start   ########################  

## scan the sample's names

sample_args <- commandArgs(trailingOnly = TRUE)
# sample_name <- sample_args[1]
# For testing: 
# sample_args <- "4C_Hb_KDM,4A_Hb_KDM,4S_Hb_KDM"
# sample_args <- "6C_Hb_KDM,6A_Hb_KDM,6S_Hb_KDM"
# sample_args <- "3C_Hb_KDM,3A_Hb_KDM,S3_Hb_KDM"   # 41 cells

sample_data <- unlist(strsplit(sample_args, ","))
sample_name_RNA <- sample_data[1]
sample_name_ATAC <- sample_data[2]
sample_name_ARC <- sample_data[3]

message('Processing sample: ', sample_name_ARC, "; RNA data in ", sample_name_RNA, " and ATAC data in ", sample_name_ATAC)


### Translate RNA barcodes - only detected cell-associated barcodes

barcode.RNA.loc_Dir <- here(cellrangerDir_GEX, sample_name_RNA, "outs", 
                            "filtered_feature_bc_matrix", "barcodes.tsv.gz")
v_sample_bc_gex <- read.csv(barcode.RNA.loc_Dir, header = FALSE)
#head(v_sample_bc_gex)
# remove GEM-Bed number `-1` and creates a vector
v_sample_bc_gex <- str_remove(unlist(v_sample_bc_gex), "-1")
df_gex <- as.data.frame(v_sample_bc_gex) # to add line-number corresponding barcode from cellranger-arc
v_sample_bc_gex <- unlist(as.data.frame(v_sample_bc_gex))

# head(v_barcodes_GEX_ARC)
# length(v_barcodes_GEX_ARC)
# head(v_sample_bc_gex)
# length(v_sample_bc_gex)

## retrieve line-number of matching vector in the GEX cellranger-arc reference
translated_bc_arc_gex <- pmatch(v_sample_bc_gex, v_barcodes_GEX_ARC)
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

translated_bc_arc_atac <- pmatch(v_sample_bc_atac, v_barcodes_ATAC_ARC)
message("Translated barcodes from cellranger-arc ATAC: ", length(translated_bc_arc_atac))

## match translated barcodes between GEX and ATAC assays alone
match_cells <- intersect(translated_bc_arc_gex, translated_bc_arc_atac)
number_match_cells <- length(match_cells)
message("Matching cells between RNA and ATAC assays: ", number_match_cells)

## Cross validate matching idx btw gex and atac
# head(sort(translated_bc_arc_gex), n=10)
# head(sort(translated_bc_arc_atac), n=10)
# head(match_cells, n=10)

## Add line-number corresponding translated ARC barcode
df_atac$ID <- translated_bc_arc_atac

## Join table with translated barcodes matching both assays
df_matching_barcodes <- inner_join(df_gex, df_atac, by = "ID") |> arrange(ID)
#df_atac[1,1]==v_barcodes_ATAC_ARC[df_atac[1,2]] & df_gex[1,1]==v_barcodes_GEX_ARC[df_gex[1,2]] #should be TRUE

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
  colnames(df_matching_barcodes) <- c("v_sample_bc_gex", "cellrangerARC_line_ID", "v_sample_bc_atac")
  write.csv(df_matching_barcodes, row.names = TRUE, quote = FALSE, 
            here(processedDir, paste0(sample_name_RNA, "_", sample_name_ATAC, "_matching_cells.csv")))
  message("Matching ARC translated barcodes saved!")

  ## Barcodes from `cellranger ARC` matching between cellranger-count and cellranger-atac 
  ## This selection of barcodes matching is for running with `cellranger-arc reanalyze` pipeline
  ARC_barcodes_matching <- v_barcodes_GEX_ARC[match_cells]
  #identical(df_matching_barcodes$v_sample_bc_gex, ARC_barcodes_matching)\ARC_barcodes_matching
  
  if (length(ARC_barcodes_matching) == sum_gex_valid_bc) {
    # write.csv(as.data.frame(ARC_barcodes_matching), row.names = FALSE, quote = FALSE, 
    #           here(cellrangerDir_reanalize, paste0(sample_name_RNA, "_matching_ARC_barcodes.csv")))
    # message("Translated barcodes for cellrangerARC reanalyze pipeline saved! (Not prefix added)")
    ## Need -1 GEM suffix to input in cellrangerARC reanalyze
    ARC_barcodes_matching_1 <- lapply(ARC_barcodes_matching, paste0, "-1")
    ARC_barcodes_matching_1 <- unlist(as.data.frame(ARC_barcodes_matching_1))
    write.csv(as.data.frame(ARC_barcodes_matching_1), row.names = FALSE, quote = FALSE, 
              here(cellrangerDir_reanalize, paste0(sample_name_RNA, "_matching_ARC_barcodes_1.csv")))
    message("Translated barcodes for cellrangerARC reanalyze pipeline saved! (with standard -1 prefix)")
    
  } else {
    message("Error saving translated barcodes for cellrangerARC reanalyze pipeline")
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

stats.data <- data.frame(
  ## samples ID
  Sample_id = paste0(sample_name_RNA, "-", sample_name_ATAC), 
  ## GEX stats
  translated_gex_bc = length(translated_bc_arc_gex), 
  valid_bc_gex_arc = sum_gex_valid_bc,
  valid_bc_gex_arc_perc = round(percent_gex, digits = 2),
  not_valid_bc_gex_arc = sum_gex_NO_valid_bc,
  not_valid_bc_gex_arc_perc = round(100-percent_gex, digits = 2),
  ## ATAC stats
  translated_atac_bc = length(translated_bc_arc_atac),
  valid_bc_atac_arc = sum_atac_valid_bc,
  valid_bc_atac_arc_perc = round(percent_atac, digits = 2),
  not_valid_bc_atac_arc = sum_atac_NO_valid_bc,
  not_valid_bc_atac_arc_perc = round(100-percent_atac, digits = 2),  
  stringsAsFactors = FALSE
)



###### Comparison of GEX_ATAC matching cells against ARC cells

## Read the raw_feature_bc_matrix.h5
message("Processing matching cells between `GEX and ATAC` against valid cell in `CellrangerARC` for sample ", sample_name_ARC) 
cellrangerDir_ARC <- here("processed-data", "cellrangerARC", sample_name_ARC, "outs", "filtered_feature_bc_matrix.h5")

## Read barcodes from ==== GEX modality =====

arc_valid_cells <- Read10X_h5(cellrangerDir_ARC) 
gex_arc <- arc_valid_cells$`Gene Expression`
v_sample_bc_ARCg <- Cells(gex_arc)
# remove GEM-Bed number `-1` and creates a vector
v_sample_bc_ARCg <- str_remove(unlist(v_sample_bc_ARCg), "-1")
#df_ARCg <- as.data.frame(v_sample_bc_ARCg) # prepare df to add line-number corresponding to cellranger-arc barcode
v_sample_bc_ARCg <- unlist(as.data.frame(v_sample_bc_ARCg))

# validation
if ( (sum((v_sample_bc_ARCg %in% v_barcodes_GEX_ARC), na.rm = TRUE) == length(v_sample_bc_ARCg)) == TRUE ) {
  translated_bc_ARCg <- pmatch(v_sample_bc_ARCg, v_barcodes_GEX_ARC)
  translated_bc_ARCg <-  translated_bc_ARCg[!is.na(translated_bc_ARCg)]
  valid_ARCg <- length(translated_bc_ARCg)
#  if ( length(v_sample_bc_ARCg) == valid_ARCg ) {
    match_cells_ARCg <- intersect(translated_bc_ARCg, match_cells)
    number_match_cells_ARCg <- length(match_cells_ARCg)
    message("Matching barcodes between cellrangerARC `GEX modality` and the matching `GEX ^ ATAC` barcodes count separately: ", number_match_cells_ARCg)
 # }
  
} else {
  message("GEX_ARC Barcodes not identified")
  valid_ARCg <- 0
  number_match_cells_ARCg <- 0
}

# ## testing complement bc in the gex sample
# length(translated_bc_arc_gex)
# length(match_cells)
# number_match_cells_ARCg
# length(intersect(translated_bc_ARCg, translated_bc_arc_gex))
# 
# v1 <-c(193, 204, 258)
# bc_arc_gex_complement <- translated_bc_arc_gex[translated_bc_arc_gex == v1] 
# bc_arc_gex_complement <- translated_bc_arc_gex[!translated_bc_arc_gex == v1] 
# match_cells_ARCg %in% translated_bc_arc_gex
# bc_arc_gex_complement <- translated_bc_arc_gex[match_cells_ARCg == translated_bc_arc_gex] 

## Read barcodes from ==== ATAC modality =====

atac_arc <- arc_valid_cells$Peaks
v_sample_bc_ARC_a <- Cells(atac_arc)
# remove GEM-Bed number `-1` and creates a vector
v_sample_bc_ARC_a <- str_remove(unlist(v_sample_bc_ARC_a), "-1")
#df_ARC_a <- as.data.frame(v_sample_bc_ARC_a) # prepare df to add line-number corresponding to cellranger-arc barcode
v_sample_bc_ARC_a <- unlist(as.data.frame(v_sample_bc_ARC_a))

# validation
# Note ATAC-ARC Barcodes are annotated with the `v_barcodes_GEX_ARC`, not with the v_barcodes_ATAC_ARC 
if ( (sum((v_sample_bc_ARC_a %in% v_barcodes_GEX_ARC), na.rm = TRUE) == length(v_sample_bc_ARC_a)) == TRUE ) {
  #translated_bc_ARC_a <- pmatch(v_sample_bc_ARC_a, v_barcodes_ATAC_ARC)
  translated_bc_ARC_a <- pmatch(v_sample_bc_ARC_a, v_barcodes_GEX_ARC)
  translated_bc_ARC_a <-  translated_bc_ARC_a[!is.na(translated_bc_ARC_a)]
  valid_ARC_a <- length(translated_bc_ARC_a)
#  if ( length(v_sample_bc_ARC_a) == valid_ARC_a) ) {
    match_cells_ARC_a <- intersect(translated_bc_ARC_a, match_cells)
    number_match_cells_ARC_a <- length(match_cells_ARC_a)
    message("Matching barcodes between cellrangerARC `ATAC modality` and the matching `GEX ^ ATAC` barcodes count separately: ", number_match_cells_ARC_a)
#  }
} else {
  message("ATAC_ARC Barcodes not identified")  
  valid_ARC_a <- 0
  number_match_cells_ARC_a <- 0
}
  
## Add stats for matching cells between cellrangerARC in both modalities and `GEX and ATAC` matching cells counted separately 
stats.data$ARC_total_cells <- length(v_sample_bc_ARCg)
#stats.data$ARC_GEX_translated_cells <- valid_ARCg
stats.data$ARC_GEX_matching_cells <- number_match_cells_ARCg
stats.data$ARC_GEX_matching_perc <- round( (number_match_cells_ARCg*100) / length(v_sample_bc_ARCg), digits=2)
#stats.data$ARC_ATAC_translated_cells <- valid_ARC_a
stats.data$ARC_ATAC_matching_cells <- number_match_cells_ARC_a
stats.data$ARC_ATAC_matching_perc <- round( (number_match_cells_ARC_a*100) / length(v_sample_bc_ARC_a), digits=2)
stats.data

write.csv(stats.data, here(processedDir, paste0(sample_name_ARC, "_", sample_name_ATAC, "_matching_cells_stats.tsv")),
            quote=FALSE, row.names = TRUE)

message(" Processed ", sample_name_RNA, " and ", sample_name_ATAC)
message("Translated barcodes from cellranger-arc ATAC: ", length(translated_bc_arc_atac))




############ Reproducibility information ####################

library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
proc.time()
options(width = 120)
session_info()

