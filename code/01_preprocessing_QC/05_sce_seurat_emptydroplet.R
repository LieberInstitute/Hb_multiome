########################################################################
## Satisfied performance of emptyDrops() in `Gene expression` multiome datasets
##
## Samples S01 t S12
## Authors. CSC + HT 
## Date. June 16th, 2023
## Last.Md Nov, 2024

## Resource:
## https://satijalab.org/seurat/archive/v3.1/conversion_vignette.html 

########################################################################

library(Seurat)
library(SingleCellExperiment)
library(DropletUtils)  
library(scales)
library(tidyverse)
library(here)
library(sessioninfo)

############        Initials        ############

here::here()

# Call functions to create and handle meta-data to seurat objects
cellrangerARC_Dir <- here("processed-data", "cellrangerARC")
processedDir_reanalyze <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze")
csvDir_reanalyze <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", "csv_files")
plotDir_reanalyze <- here("plots", "01_preprocessing_QC", "cellrangerARC_reanalyze")

# Check directories exists
if (!dir.exists(csvDir_reanalyze)) { dir.create(csvDir_reanalyze) }

#source(here("code", "remote_seurat_functions.R"))

# Define theme function
define_theme <- function(size = 15) {
  theme_bw() +
    theme(text = element_text(size = size))
}

# Load raw barcode data or Seurat with raw data
b_h5file <- TRUE

## commandArgs scans the arguments which have been supplied when the current R script was invoked (from shell sh)
sample_tmp <- commandArgs(trailingOnly = TRUE)
# For testing: 
# sample_tmp <- "S3_Hb_KDM_reanalysis, S3_Hb_KDM"

sample_data = unlist(strsplit(sample_tmp,","))
sample_name <- trimws(sample_data[[2]])

message("Reading CellRangerARC sample: ", sample_name)


# Specify if raw data will be load directly from a H5 file or a SeuratOBJ
if (!b_h5file) {

     # # The raw file is selected instead of the truncared filtered matrix
     # s_featured_bc_mtx <- here("raw-data/PBMC_CellSorted_ARC2_0_0", "pbmc_granulocyte_sorted_10k_raw_feature_bc_matrix.h5")
     # #s_meta_fname <- here("raw-data/PBMC_CellSorted_ARC2_0_0", "pbmc_granulocyte_sorted_3k_per_barcode_metrics.csv")
     # SeuratOBJ <- f_create_seurat_RAW(s_sample_name2, s_featured_bc_mtx)
     # #s_frag_namefile <- here('raw-data/PBMC_CellSorted_ARC2_0_0/', 'pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz')
     # # Make a Seurat obj a SCE
     # DefaultAssay(SeuratOBJ) <- "RNA"
     # SeuratOBJ.sce <- as.SingleCellExperiment(SeuratOBJ)

 } else if (b_h5file) {

    # Load h5 file directly as usually do in single cell exp
    cellrangerARC_Dir <- here(cellrangerARC_Dir,  sample_name, "outs", "raw_feature_bc_matrix.h5")  
    h5_raw_path <- Read10X_h5(here(cellrangerARC_Dir)) 
    raw.sce <- h5_raw_path$`Gene Expression`

} 

#str(raw.sce)
head(raw.sce, n=3)

totalCells <- length(Cells(raw.sce))

message("Processing ", totalCells, " from CellRanger ARC sample ", sample_name)
# Processing 583052 from multiome sample 4S_Hb_KDM

# barcodeRanks method (from DropletUtils package), compute barcode rank statistics and identify the knee and inflection points on the total count curve
#   @fit.bounds specify the lower and upper bounds on the total UMI count from which to obtain a section of the curve for spline fitting
bcRanks <- barcodeRanks(raw.sce, fit.bounds = c(10, 1e3))

# plot the computed barcode rank. The rank of each barcode (averaged across ties).
range((bcRanks$rank))   # [1]      1.0 441,982.5

range((bcRanks$total))  # [1]     0 93,817

knee_lower <- metadata(bcRanks)$knee + 100
# [1] 998

infection <- metadata(bcRanks)$inflection
# [1] 276

o <- order(bcRanks$rank)
# legend_label <- c(paste0("knee: ", toString(knee_lower)), paste0("inflection: ", toString(infection)))
# plot(bcRanks$rank, bcRanks$total, log="xy", xlab="Rank", ylab="Total UMI Counts") +
#   lines(bcRanks$rank[o], bcRanks$fitted[o], col="red") +
#   abline(h=metadata(bcRanks)$knee, col="dodgerblue", lty=2) +
#   abline(h=metadata(bcRanks)$inflection, col="forestgreen", lty=2) +
#   title(paste0(sample_name, '\nWithout gene filtering')) +
#   legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
#        legend = legend_label)

#### Run emptyDrops w/ knee + 100 ####
set.seed(05112024)

# emptyDrops from DropletUtils Bioconductor package is used
# Distinguish between droplets containing cells and ambient RNA in a droplet-based single-cell RNA sequencing experiment. 
sce.out <- DropletUtils::emptyDrops(
  raw.sce,
  niters = 30000,       # number of iterations to use for the Monte Carlo p-value calculations
  lower = knee_lower    # numeric scalar specifying the lower bound on the total UMI count, below this all barcodes are assumed empty droplets
)

# Saving data
save(sce.out, file = here(processedDir_reanalyze, paste0(sample_name, "_droplet_scores.rds")))
# ~/processed-data/01_preprocessing_QC/cellrangerARC_reanalyze/
message("Saved droplet scores RDS object")

#### QC Plot ####

FDR_cutoff <- 0.001
table(Signif = sce.out$FDR <= FDR_cutoff)
# Signif
# FALSE  TRUE 
# 311  7241 

#  add arbitrary margins on a multidimensional array
tab <- addmargins(table(Signif = sce.out$FDR <= FDR_cutoff, 
                 Limited = sce.out$Limited, 
                 useNA = "ifany"))
write.csv(tab, here(csvDir_reanalyze, paste0(sample_name, "_droplet_FDRcutoff.csv")))
n_cell_anno <- paste("Non-empty:", sum(sce.out$FDR < FDR_cutoff, na.rm = TRUE))

# Calculate non Emptydroplets value
nonEmptydroplets <- (sce.out |> as.data.frame() |> filter(FDR < FDR_cutoff) |> summarise(n = n()))$n

# Calculate percentage of non Emptydroplets
per.nonemptydroplets <- ((nonEmptydroplets*100) / totalCells) #, digits = 4)

# create results dataframe
# https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/03_droplet_qc.R

data = data.frame(Cells=c("TotalCells", "NonEmptyCells", "PercentageNonEmpty"), 
                  values=c(totalCells, nonEmptydroplets, paste(per.nonemptydroplets, "%")))
message(paste0('Total True Cells on Sample ', sample_name))
print(data)
write.csv(data, here(csvDir_reanalyze, paste0(sample_name, "_total_TrueCells.csv")))
# Cells             values
# 1         TotalCells             583052
# 2      NonEmptyCells               7241
# 3 PercentageNonEmpty 1.24191324272964 %

head(sce.out, n=3)
#rowsum(sce.out$FDR>0)
# <integer> <numeric> <numeric> <logical> <numeric>
# AAACAGCCAAACAACA-1         2        NA        NA        NA        NA
# AAACAGCCAAACATAG-1         1        NA        NA        NA        NA
# AAACAGCCAAACCCTA-1         1        NA        NA        NA        NA


# Prepare data frame with additional FDR column

droplet_elbow_data <- as.data.frame(bcRanks) |> mutate(FDR = sce.out$FDR)

## Define parameters to plot
knee_meta <- metadata(bcRanks)$knee
knee_lower_label <- paste0("Knee est 'lower' (", knee_lower, ')')
second_knee_label <- paste0("Second Knee (", knee_meta, ')') 
title <- paste0("Sample: ", sample_name)
subtitle <- n_cell_anno

## Create ggplot object
message("Creating elbow plot ...")

droplet_elbow_plot <- droplet_elbow_data |>
  ggplot(aes(x = rank, y = total, color = FDR < FDR_cutoff)) +
  # Define points
  geom_point(alpha = 0.5, size = 1) +
  # Define lines and annotations
  geom_hline(yintercept = knee_meta, linetype = "dotted", color = "gray") +
  annotate("text", x = 10, y = knee_meta, label = second_knee_label, vjust = -1, color = "gray") +
  geom_hline(yintercept = knee_lower, linetype = "dashed") +
  annotate("text", x = 10, y = knee_lower, label = knee_lower_label, vjust = -0.5) +
  # Define scales
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10", oob = scales::squish_infinite) +
  
  # Define labels
  labs(
    x = "Barcode Rank",
    y = "Total UMI Counts",
    title = title,
    subtitle = subtitle,
    color = paste("FDR <", FDR_cutoff)
  ) +
  # Apply theme
  define_theme() +
  theme(legend.position = "bottom")

## Save the png format
plt_name <- here(plotDir_reanalyze, paste0(sample_name, "_droplet_qc.png"))
ggsave(filename = plt_name, plot = droplet_elbow_plot)

message("Saved elbow plot!")

message("Done!")


# library(slurmjobs)
# job_single(
#   name = "05_sce_seurat_emptydroplet", memory = "30G", cores = 1, create_shell = TRUE,
#   task_num = 10
# )


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
