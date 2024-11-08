########################################################################
## Applying QC metrics to Seurat object
##
## Note. Compute the QC metrics with scran and scater 
########################################################################

library("SingleCellExperiment")
library("scuttle")
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

# load(here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", paste0(Seurat_base_name, "_droplet_scores.rds"))) # sce.out

## Prepate dirs and read the raw matrix 
cellrangerDir_reanalyze <- here("processed-data", "cellrangerARC", Seurat_base_name, "outs")  # to process integrated Seurat Object 
processedDir_reanalyze <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze")
plotDir_reanalyze <- here("plots", "01_preprocessing_QC", "cellrangerARC_reanalyze")

## Load raw data + droplets results
unfiltered_path <- here(cellrangerDir_reanalyze, "raw_feature_bc_matrix.h5")
sce_emptydrops_path <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", paste0(Seurat_base_name, "_droplet_scores.rds")) # sce.out

# Check processed_data and plot directories exists
if (!dir.exists(processedDir_reanalyze)) { dir.create(processedDir_reanalyze) }
if (!dir.exists(plotDir_reanalyze)) { dir.create(plotDir_reanalyze) }

## Read the raw_feature_bc_matrix.h5
message("Reading raw feature bc data corresponding to sample: ", Seurat_base_name) # ../cellrangerARC/S1_Hb_KDM/outs/raw_feature_bc_matrix.h5
# sce.raw <- read10xCounts(here(cellrangerDir_reanalyze, "raw_feature_bc_matrix.h5")) # dgCMatrix data. Barcodes for columns and genes by rows (DropletUtils)
# sce.raw

set.seed(5112024)


# ############### FINDING HIGH MITO ##############################################

## Call function
process_sample <- function(sample_path, sce_out_path, fdr_threshold = 0.001) {
  # fdr_threshold = 0.001
  # raw_sample_path <- here(cellrangerDir_reanalyze, "raw_feature_bc_matrix.h5")
  # sce_out_path <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", paste0(Seurat_base_name, "_droplet_scores.rds"))
  
  # Load the sample data
  sce <- read10xCounts(raw_sample_path)
  
  # Unifying feature names
  rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
  
  # Map the IDs
  location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, column="SEQNAME", keytype="GENEID")
  
  # Load sce.out
  load(sce_out_path, verbose = T)
  
  # Subset our SingleCellExperiment object to retain only the detected cells
  sce <- sce[,which(sce.out$FDR <= fdr_threshold)]
  unfiltered <- sce
  total_unfiltered_cells <- dim(assay(unfiltered))[2]
  
  # Quality control
  # Filtering on the mitochondrial proportion
  is.mito <- grep("MT-", rownames(sce))
  stats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
  # colnames(stats)
  # summary(stats$subsets_Mito_percent)
  # stats <- perCellQCMetrics(sce, subsets=list(Mito=which(location=="MT")))
  high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
  sce <- sce[,!high.mito]
  total_filtered_cells <- dim(assay(sce))[2]
  #summary(high.mito)
  
  ## store this in the colData() of our SingleCellExperiment object for future reference
  colData(unfiltered) <- cbind(colData(unfiltered), stats)
  unfiltered$discard <- high.mito
  # colnames(colData(unfiltered))

  ## Build title labels
  out_detected <- total_unfiltered_cells - total_filtered_cells
  out_detected_p <- round( ((out_detected*100) / total_unfiltered_cells), digits = 2 ) 
  caption_label <- paste0(out_detected, " cells (", out_detected_p, "%) outliers detected from ", total_unfiltered_cells, ". ", total_filtered_cells, " True cells.")
  ## Build plot
  plot_grid <- gridExtra::grid.arrange(
    plotColData(unfiltered, y="sum", colour_by="discard") +
      scale_y_log10() + ggtitle("Total count"),
    plotColData(unfiltered, y="detected", colour_by="discard") +
      scale_y_log10() + ggtitle("Detected features"),
    plotColData(unfiltered, y="subsets_Mito_percent",
                colour_by="discard") + ggtitle("Mito percent"),
    ncol = 3,
    top = paste0(Seurat_base_name, " Outliers detected"),
    bottom = caption_label
  )

  # Save the plot
  plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_isOutliers_metrics.png"))
  #ggsave(filename = plotName, plot = plot_grid)
  
  message("Done!")
  
  return(sce)

}

## call function to process sample 
process_sample(unfiltered_path, sce_emptydrops_path)


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
