########################################################################
## Applying QC metrics to filtered Seurat object
##
## Note. Compute the QC metrics with Seurat, Scran and Scater
########################################################################


# Load libraries
library("Seurat")
library("SingleCellExperiment")
#library("VariantAnnotation")
library("here")
library("ggplot2")
library("ggrepel")
library("scater")
library("batchelor")
library("scran")
library("scry")
library("uwot")
#library("DropletUtils")
#library("Rtsne")
library("gridExtra")
#library("EnsDb.Hsapiens.v86")
#library("reshape")
#library("cowplot")
library("dplyr")
library("sessioninfo")


## Read directories

here::here()

#cellrangerDir_reanalyze <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze") # to process individual Seurat Objects
cellrangerDir_reanalyze <- here("processed-data", "02_merge_seurats", "cellrangerARC_reanalyze")  # to process intgerated Seurat Object 
processedDir_reanalyze <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze")
plotDir_reanalyze <- here("plots", "01_preprocessing_QC", "cellrangerARC_reanalyze")

# Check processed_data and plot directories exists
if (!dir.exists(processedDir_reanalyze)) { dir.create(processedDir_reanalyze) }
if (!dir.exists(plotDir_reanalyze)) { dir.create(plotDir_reanalyze) }

## function to load pre-existing Seurat
get_seurat <- function(name) { sobj <- readRDS(name); return(sobj)}

## Build Seurat object name
Seurat_base_name <- 'seurat.norm_counts' 
rds_name <- here(cellrangerDir_reanalyze, paste0(Seurat_base_name, '_Harmony.rds'))
#rds_name

## Load Seurat object
SeuratOBJ <- get_seurat(rds_name)

## Verify Seurat object
table(SeuratOBJ$orig.ident)

message("Processing QCs for ", Seurat_base_name)





set.seed(777)

# Loading droplets results 
load(here("processed-data", "08_build_sce", "droplet_scores_Hippo_42_4.RDS"), verbose = T)
sce.out
# QC approach based on Workflow 3.3 in OSCA:
# http://bioconductor.org/books/3.16/OSCA.workflows/unfiltered-human-pbmcs-10x-genomics.html#quality-control-2

############### FINDING HIGH MITO ##############################################
# read the raw data matrix
sce.42_4 <- read10xCounts(here("rafael_rerun/42_4/outs", "raw_feature_bc_matrix.h5"))
# class: SingleCellExperiment 
# dim: 128902 730781 
# metadata(1): Samples
# assays(1): counts
# rownames(128902): ENSG00000243485 ENSG00000237613 ...
# KI270713.1:31342-32216 KI270713.1:34152-35040
# rowData names(3): ID Symbol Type
# colnames: NULL
# colData names(2): Sample Barcode
# reducedDimNames(0):
#     mainExpName: NULL
# altExpNames(0):
#     
# Unifying feature names
rownames(sce.42_4) <- uniquifyFeatureNames(
  rowData(sce.42_4)$ID, rowData(sce.42_4)$Symbol)

location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.42_4)$ID, 
                   column="SEQNAME", keytype="GENEID")
head(location, n=10)
# Warning message:
#   Unable to map 94912 of 128902 requested IDs. 

# Once we are satisfied with the performance of emptyDrops(), we subset our SingleCellExperiment object to retain only the detected cells. Discerning readers will notice the use of which(), which conveniently removes the NAs prior to the subsetting

sce.42_4 <- sce.42_4[,which(sce.out$FDR <= 0.001)]
unfiltered <- sce.42_4

# Quality control
# Filtering on the mitochondrial proportion
stats <- perCellQCMetrics(sce.42_4, subsets=list(Mito=which(location=="MT")))

# Setup parameters to state different levels of outliers 
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
sce.42_4 <- sce.42_4[,!high.mito]
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
  sum = sce.42_4.qc@listData$sum,
  subsets_MT_percent = sce.42_4.qc@listData$subsets_MT_percent
)

# Defining thresholds
thresholds <- attr(discard.mito, "thresholds")["higher"]

# Creating the plot
ggplot(data = data, aes(x = sum, y = subsets_MT_percent)) +
  geom_point() +
  scale_x_log10() +
  geom_hline(yintercept = thresholds, color = "red") +
  labs(
    title = "Scatter plot of Total count vs. Mitochondrial %",
    x = "Total count",
    y = "Mitochondrial %"
  ) +
  theme_minimal()
