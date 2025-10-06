########################################################################
## Summaries/Prepare plots for interpretate Peaks Overlapping Links-DARs
##
## Authors. CSC
## Date. Oct 06, 2025
## Recommended resources on interactive mode: srun --pty --mem=15GB --x11 bash
########################################################################

library("dplyr")
library("purrr")
library("ggplot2")
library("patchwork")
library("tidyverse")
library("tidyr")
library("stringr")
library("here")

#===============================================================================
# resolution_level = "Mid"      # 18 cell-types
#===============================================================================

## setup variable names

# resolution_level = "Mid" 
FDR = 0.2 # actual value to run the DAR-Links
# we do not use log FC for this exploratory analysis

## Set directory names
inputCSV_Overlaps_Dir <- here(
    "processed-data",
    "06_peak_calling",
    "20_explore_overlapping_linkpeaks_DARs"
)
processedDir <- here(
    "processed-data",
    "06_peak_calling",
    "21_explore_top_overlapping_linkpeaks_DARs"
)
plotDir <- here(
    "plots",
    "06_peak_calling",
    "21_explore_top_overlapping_linkpeaks_DARs"
)

if (!dir.exists(processedDir)) { dir.create(processedDir) }
if (!dir.exists(plotDir)) { dir.create(plotDir) }

##==============================================================================

message("Loading Unique-Overlap with directionality ...")

## load full overlaps df
overlaps_df <- read.csv(here(inputCSV_Overlaps_Dir, 
                             paste0("overlaps_unique_with_direction_ct_FDR", FDR, ".csv")))
## test
colnames(overlaps_df)
nrow(overlaps_df) # 3201




message("All plots done!!!")


# library("slurmjobs")
# job_single(
#   "19_Linkage_DARs_analysis",
#   create_shell = TRUE,
#   partition = "katun",
#   memory = "30G",
#   cores = 2,
#   logdir = "logs",
#   command = "Rscript 19_Linkage_DARs_analysis.R",
#   create_logdir = FALSE
# )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



