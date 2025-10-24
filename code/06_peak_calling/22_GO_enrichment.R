########################################################################
## GO enrichment in cCRE and classified Linked OCR and DARs unlinked
##
## Authors. CSC
## Date. Oct 24, 2025
## Recommended resources on interactive mode: srun --pty --mem=20GB --x11 bash
########################################################################

#library("pheatmap")
#library("reshape2")
library("dplyr")
library("purrr")
library("ggplot2")
library("patchwork")
#library("ggrepel")
library("tidyverse")
library("tidyr")
library("stringr")
library("here")

library("getopt")
library("org.Hs.eg.db")
library("clusterProfiler")
library("rrvgo")
library("ComplexHeatmap")

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================


## Set directory names
inputCSV_cCRE_ORC_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "21_overlaping_FDRscores_TopHeatmap"
)
processedDir <- here(
    "processed-data",
    "06_peak_calling",
    "22_GO_enrichment"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "22_GO_enrichment"
)

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }


##==============================================================================

message("Loading cCRE and ORC file ...")

f_name = "overlaps_linkPeak_DARs_classified_thr_CC0.3_thr_DAR0.1.csv"

## load raw overlaps and 2-shared ct overlaps 
cCRE_df <- read.csv(here(inputCSV_cCRE_ORC_Dir, f_name))
nrow(cCRE_df) # 5603
head(cCRE_df)

## validation
unique(cCRE_df$type_classification)
as.data.frame(table(cCRE_df$category))

only_cCRE_df <- cCRE_df |> 
    filter(category %in% c("cell-specific cCRE (-)", "cell-specific cCRE (+)"))

table(only_cCRE_df$category)
# cell-specific cCRE (-) cell-specific cCRE (+) 
# 657                    706 
table(only_cCRE_df$overlap_type)
# Shared Unique 
# 508    855 

cluster_levels <- only_cCRE_df$cell_type |> unique()






