########################################################################
## Applying QC metrics to filtered Seurat object
##
## Note. Compute the QC metrics with Seurat, Scran and Scater
########################################################################

# Load libraries
#library("Seurat")
#library("Signac")
#library("scuttle")
#library("VariantAnnotation")
#library("Rtsne")
#library("reshape")
#library("cowplot")
#library("dplyr")

library("SingleCellExperiment")
library("here")
library("ggplot2")
library("ggrepel")
library("scater")
library("batchelor")
library("scran")
library("scry")
library("uwot")
library("DropletUtils")
library("gridExtra")
library("EnsDb.Hsapiens.v86")
library("sessioninfo")

## Read directories

here::here()

# test
Seurat_base_name <- "4S_Hb_KDM"

## Loading droplets results 
load(here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", paste0(Seurat_base_name, "_droplet_scores.rds"))) 
## Read the raw matrix 
cellrangerDir_reanalyze <- here("processed-data", "cellrangerARC", Seurat_base_name, "outs")  # to process integrated Seurat Object 
processedDir_reanalyze <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze")
plotDir_reanalyze <- here("plots", "01_preprocessing_QC", "cellrangerARC_reanalyze")

# Call functions to create and handle meta-data to seurat objects
#source(here("code/functions_custom", "remote_file_caller.R"))

# Check processed_data and plot directories exists
if (!dir.exists(processedDir_reanalyze)) { dir.create(processedDir_reanalyze) }
if (!dir.exists(plotDir_reanalyze)) { dir.create(plotDir_reanalyze) }

## Read the raw_feature_bc_matrix.h5
message("Reading raw feature bc data for sample: ", Seurat_base_name) # ../cellrangerARC/S1_Hb_KDM/outs/raw_feature_bc_matrix.h5
sce.raw <- read10xCounts(here(cellrangerDir_reanalyze, "raw_feature_bc_matrix.h5")) # dgCMatrix data. Barcodes for columns and genes by rows (DropletUtils)
sce.raw

# # Extract the 'Gene Expression' matrix only
# raw.sce <- h5_raw_path$`Gene Expression` # SingleCellExperiment data.
# 
# ## Get total number of raw cells in the gene expression assay
# totalCells <- length(Cells(raw.sce))
# # [1] 583052
# 
# message("Building Seurat object with raw data for ", Seurat_base_name)
# 
# # ## Conversion from SingleCellExperiment objects to Seurat objects
# # typeof(SeuratOBJ) 
# # SeuratOBJ@assays$RNA@layers$counts
# # #sce <- SingleCellExperiment(list(counts=as.matrix(SeuratOBJ@assays$RNA@layers$counts))) # Only pull counts
# # seurat.sce <- as.SingleCellExperiment(SeuratOBJ) # pull all assays in both RNA and ATAC assays
# # dim(seurat.sce)
# # # [1] 70990  7345
# # colnames(colData(seurat.sce))

set.seed(777)

#### Compute QC metrics ####

# # Initialize the Seurat object with the raw (non-normalized data)
# Seurat.raw <- CreateSeuratObject(
#   counts = raw.sce,
#   assay = "RNA",
#   project = Seurat_base_name
#   #meta.data = meta2
# )
# Seurat.raw
# # An object of class Seurat 
# # 36601 features across 583052 samples within 1 assay 
# # Active assay: RNA (36601 features, 0 variable features)
# # 1 layer present: counts
# str(Seurat.raw)
# Seurat.sce <- SingleCellExperiment(list(counts=as.matrix(Seurat.raw@assays$RNA@layers$counts))) # Only pull counts
# str(Seurat.sce)
#  
# 
# ## Read filtered feature bc 
# sample_path <- paste0(here(processedDir_reanalyze, Seurat_base_name), "_reanalysis.rds")
# Seurat.filter <- readRDS(sample_path)
# Seurat.filter
# str(Seurat.filter)
# Seurat.sce.f <- SingleCellExperiment(list(counts=as.matrix(Seurat.filter@assays$RNA@layers$counts)))
# str(Seurat.sce.f)
# 
# # Then we can add an simple CPM transformation to the original matrix count matrix and store it
# exprs(Seurat.sce.f) <- log2(calculateCPM(Seurat.sce.f) + 1)  #SCATER



# ############### FINDING HIGH MITO ##############################################
# # read the raw data matrix
# sce.42_4 <- read10xCounts(here("rafael_rerun/42_4/outs", "raw_feature_bc_matrix.h5"))

# Unifying feature names
rownames(sce.raw) <- uniquifyFeatureNames(rowData(sce.raw)$ID, rowData(sce.raw)$Symbol)

location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.raw)$ID, column="SEQNAME", keytype="GENEID")
head(location, n=10)
# Warning message:
#   Unable to map 73631 of 107621 requested IDs.

# We subset our SingleCellExperiment object to retain only the detected cells. Discerning readers will notice the use of which(), which conveniently removes the NAs prior to the sub-setting

sce.raw <- sce.raw[,which(sce.raw$FDR <= 0.001)]
unfiltered <- sce.raw

# Quality control
# Filtering on the mitochondrial proportion
stats <- perCellQCMetrics(sce.raw, subsets=list(Mito=which(location=="MT")))

# Setup parameters to state different levels of outliers 
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
sce.raw <- sce.raw[,!high.mito]
summary(high.mito)
# Mode   FALSE    TRUE 
# logical    8429    1441 


colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- high.mito

gridExtra::grid.arrange(
  plotColData(unfiltered, y="sum", colour_by="discard") +
    scale_y_log10() + ggtitle("Total count"),
  plotColData(unfiltered, y="detected", colour_by="discard") +
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(unfiltered, y="subsets_Mito_percent",
              colour_by="discard") + ggtitle("Mito percent"),
  ncol=2
)

# Write a function
process_sample <- function(sample_path, sce_out_path, output_path, fdr_threshold = 0.001) {
  
  # Load the sample data
  sce <- read10xCounts(sample_path)
  
  # Unifying feature names
  rownames(sce) <- uniquifyFeatureNames(
    rowData(sce)$ID, rowData(sce)$Symbol)
  
  # Map the IDs
  location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, 
                     column="SEQNAME", keytype="GENEID")
  
  # Load sce.out
  load(sce_out_path, verbose = T)
  
  # Subset our SingleCellExperiment object to retain only the detected cells
  sce <- sce[,which(sce.out$FDR <= fdr_threshold)]
  unfiltered <- sce
  
  # Quality control
  # Filtering on the mitochondrial proportion
  stats <- perCellQCMetrics(sce, subsets=list(Mito=which(location=="MT")))
  high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
  sce <- sce[,!high.mito]
  summary(high.mito)
  
  colData(unfiltered) <- cbind(colData(unfiltered), stats)
  unfiltered$discard <- high.mito
  
  plot_grid <- gridExtra::grid.arrange(
    plotColData(unfiltered, y="sum", colour_by="discard") +
      scale_y_log10() + ggtitle("Total count"),
    plotColData(unfiltered, y="detected", colour_by="discard") +
      scale_y_log10() + ggtitle("Detected features"),
    plotColData(unfiltered, y="subsets_Mito_percent",
                colour_by="discard") + ggtitle("Mito percent"),
    ncol=2
  )
  
  # Save the plot
  ggsave(filename = output_path, plot = plot_grid)
  
  return(sce)
}

# Define the paths
# 42_4
sample_42_4_path = here("rafael_rerun/42_4/outs", "raw_feature_bc_matrix.h5")
sce_out_42_4_path = here("processed-data", "08_build_sce", "droplet_scores_Hippo_42_4.RDS")
# 42_1
sample_42_1_path = here("rafael_rerun/42_1/outs", "raw_feature_bc_matrix.h5")
sce_out_42_1_path = here("processed-data", "08_build_sce", "droplet_scores_Hippo_42_1.RDS")

# Define the plots path
output_42_4_path = here("plots/09_sce_qcmetter", "42_4_qc_metrics.png")
output_42_1_path = here("plots/09_sce_qcmetter", "42_1_qc_metrics.png")

# Run the function
sce_42_4 <- process_sample(sample_42_4_path, sce_out_42_4_path, output_42_4_path)
sce_42_1 <- process_sample(sample_42_1_path, sce_out_42_1_path, output_42_1_path)

# Creating a new data frame from the given S4 object lists

data <- data.frame(
  #sum = sce.42_4.qc@listData$sum,
  sum = sce.out@listData$sum,
  #subsets_MT_percent = sce.42_4.qc@listData$subsets_MT_percent
  subsets_MT_percent = sce.out@listData$subsets_MT_percent
)

# Defining thresholds
thresholds <- attr(discard.mito, "thresholds")["higher"]

# Creating the plot
ggplot(data = data, aes(x = sum, y = subsets_MT_percent)) +
  geom_point() +
  scale_x_log10() +
  #geom_hline(yintercept = thresholds, color = "red") +
  labs(
    title = "Scatter plot of Total count vs. Mitochondrial %",
    x = "Total count",
    y = "Mitochondrial %"
  ) +
  theme_minimal()
