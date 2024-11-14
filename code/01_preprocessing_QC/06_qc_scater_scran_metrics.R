########################################################################
## Calculate standard quality control metrics based on the singleCellExperiment assay in Seurat's GEX assays
## CSC. Nov-2024
########################################################################

library("SingleCellExperiment")
library("Seurat")
library("scuttle")
library("here")
library("ggplot2")
# library("ggrepel") # provides geoms for ggplot2 to repel overlapping text labels
library("scater")
# library("batchelor") # Implements a variety of methods for batch correction of single-cell (RNA sequencing) data
library("scran")
library("scry")
# library("uwot") # An R implementation of the Uniform Manifold Approximation and Projection (UMAP)
library("DropletUtils")
library("gridExtra")
library("EnsDb.Hsapiens.v86")
library("sessioninfo")

## Read directories

here::here()

## Scans arguments invoked from slurm job shell sh
sample_tmp <- commandArgs(trailingOnly = TRUE)
# For testing:
# sample_tmp <- "4S_Hb_KDM_reanalysis, 4S_Hb_KDM"
sample_data = unlist(strsplit(sample_tmp,","))
Seurat_base_name <- trimws(sample_data[[2]])

message("Reading CellRangerARC reanalyze sample: ", Seurat_base_name)

## Prepare Dir(s) and read the raw matrix 
cellrangerDir_reanalyze <- here("processed-data", "cellrangerARC", Seurat_base_name, "outs")  
csvDir_reanalyze <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", "csv_files")
plotDir_reanalyze <- here("plots", "01_preprocessing_QC", "cellrangerARC_reanalyze")

## Load raw data
unfiltered_path <- here(cellrangerDir_reanalyze, "raw_feature_bc_matrix.h5")
# sce_emptydrops_path <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", paste0(Seurat_base_name, "_droplet_scores.rds")) # sce.out

# Check processed_data and plot directories exists
if (!dir.exists(csvDir_reanalyze)) { dir.create(csvDir_reanalyze) }
if (!dir.exists(plotDir_reanalyze)) { dir.create(plotDir_reanalyze) }

## Read the raw_feature_bc_matrix.h5
message("Reading raw feature barcode data corresponding to sample: ", Seurat_base_name) # ../cellrangerARC/S1_Hb_KDM/outs/raw_feature_bc_matrix.h5

set.seed(11112024)


# ############### FINDING HIGH MITO ##############################################

## Call function
process_sample <- function(sample_path, sce_out_path, fdr_threshold = 0.001) {
  # fdr_threshold = 0.001   # Only need it if you are processing empyDrops results
  # sample_path <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/cellrangerARC/4S_Hb_KDM/outs/raw_feature_bc_matrix.h5"
  
  # Load the sample data
  sce <- read10xCounts(sample_path)
  # class: SingleCellExperiment 
  # dim: 181201 702703 
  
  # Unifying feature names
  rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
  # Map the IDs
  location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, column="SEQNAME", keytype="GENEID")
  
  ## If interested in use the emptyDrops derived cells for calculate Outliers use the chunk code below
  # sce_out_path <- here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", paste0(Seurat_base_name, "_droplet_scores.rds"))
  # load(sce_out_path, verbose = T) # load sce.out (emptyDrops derived RDS object)
  ## Subset our SingleCellExperiment object to retain only the detected cells
  # sce <- sce[,which(sce.out$FDR <= fdr_threshold)]
  # unfiltered <- sce
  
  ## Load sce.out from cellrangerARC_reanalyze, then transform to sce object
  sce.out2 <- readRDS(here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", paste0(Seurat_base_name, "_reanalysis.rds")))
  sce.out2 <- as.SingleCellExperiment(sce.out2, assay = "RNA")

  ## Select those cells detected from custom Cell RangerARC reanalyze pipeline
  v_reanalyze_cells <- Cells(sce.out2)

  ## Subset to keep only valid barcodes from Cell RangerARC reanalyze
  sce <- sce[,sce$Barcode %in% c(v_reanalyze_cells)]
  # dim(sce) #[1] 181201   5050
  #unfiltered <- sce
  
  total_unfiltered_cells <- ncol(sce) # cells in cols
  
  ## Quality control, check low quality cells
  
  #sce <- scuttle::logNormCounts(sce) ?
  
  ## High mito
  is.mito <- grep("MT-", rownames(sce))
  sce <- scuttle::addPerCellQC(
    sce,
    subsets = list(Mito = is.mito),
    BPPARAM = BiocParallel::MulticoreParam(4)
  )
  
  ## For reference, these are result derived from reanalyze (5050 cells) vs emptyDrops for sample S03 (5802 cells)
  #summary(sce$subsets_Mito_percent) 
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.0000  0.4765  1.0858  2.5113  2.1428 84.0290
  
  sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher") # batch = sce$Sample, running one sample at time
  ## cells pass
  table(sce$high_mito)
  ## store the stats in the colData() for future reference
  # colData(unfiltered) <- cbind(colData(unfiltered), stats)
  # unfiltered$discard <- high.mito
  # colnames(colData(unfiltered))

  ## low library size
  sce$low_sum <- isOutlier(sce$sum, log = TRUE, type = "lower") # , batch = sce$Sample
  table(sce$low_sum)
  # FALSE 
  # 5050 
  
  ## low detected features
  sce$low_detected <- isOutlier(sce$detected, log = TRUE, type = "lower") #, batch = sce$Sample
  table(sce$low_detected)
  # FALSE 
  # 5050 
  
  ## All low sum are also low detected
  table(sce$low_sum, sce$low_detected)
  
  ## Annotate cells to remove
  sce$discard_auto <- sce$high_mito | sce$low_sum | sce$low_detected
  
  table(sce$discard_auto)
  # FALSE  TRUE 
  # 4523   527 
  
  # ## discard 10.43 % for S03
  # 100 * sum(sce$discard_auto) / ncol(sce)
  # (qc_t <- addmargins(table(sce$Sample, sce$discard_auto)))
  # ## in percentage
  # round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)

  # Filter cells that PASS Outliers and save barcodes filtered
  sce_bc <- sce[,!sce$discard_auto]
  ncol(sce_bc)
  total_filtered_cells <-  length(sce_bc$discard_auto[sce$discard_auto==FALSE])
  
  message("Total cells filtered (PASS) from sample ", Seurat_base_name, ": ", total_filtered_cells, " from ", total_unfiltered_cells)

  csv_name <- here(csvDir_reanalyze, paste0(Seurat_base_name, "_bc_PASS_isOutliers.csv"))
  write.csv(sce_bc$Barcode, csv_name)
  
  message(paste0("Saved valid (PASS) barcodes for sample ", Seurat_base_name))
  
  ## Build title labels
  out_detected <- total_unfiltered_cells - total_filtered_cells
  out_detected_p <- round( ((out_detected*100) / total_unfiltered_cells), digits = 2 ) 
  caption_label <- paste0(out_detected, " cells (", out_detected_p, "%) with outliers detected. Filtered ", total_filtered_cells, " from ", total_unfiltered_cells)

  ## Build plot
  plot_grid <- gridExtra::grid.arrange(
    plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "high_mito") + ggtitle("Mito Precent"), 
    ## low sum
    plotColData(sce, x = "Sample", y = "sum", colour_by = "low_sum") + scale_y_log10() + ggtitle("Total count"),    
    ## low genes
    plotColData(sce, x = "Sample", y = "detected", colour_by = "low_detected") + scale_y_log10() + ggtitle("Detected features"), 
    ncol = 3,
    top = paste0("Outliers detected. Sample ", Seurat_base_name),
    bottom = caption_label
  )

  # Save the plot
  plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_isOutliers_metrics.png"))
  ggsave(filename = plotName, plot = plot_grid)
  print(paste0("Saved plot from isOutliers for sample ", Seurat_base_name))
  
  # # Mito rate vs n detected features
  # plotColData(sce,
  #             x = "detected", y = "subsets_Mito_percent",
  #             colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5)
  # # Detected features vs total count
  # plotColData(sce,
  #             x = "sum", y = "detected",
  #             colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5)
  
  
  return(sce)

}

## call function to process sample 
process_sample(unfiltered_path, sce_emptydrops_path)

message("Done!")


# # Creating a new data frame from the given S4 object lists
# 
# sce.out@listData$sum
# 
# data <- data.frame(
#   sum = sce.out@listData$sum,
#   subsets_MT_percent = sce.out@listData$subsets_MT_percent
# )
# 
# # Defining thresholds
# thresholds <- attr(discard.mito, "thresholds")["higher"]
# 
# # Creating the plot
# ggplot(data = data, aes(x = sum, y = subsets_MT_percent)) +
#   geom_point() +
#   scale_x_log10() +
#   #geom_hline(yintercept = thresholds, color = "red") +
#   labs(
#     title = "Scatter plot of Total count vs. Mitochondrial %",
#     x = "Total count",
#     y = "Mitochondrial %"
#   ) +
#   theme_minimal()


## Add job array
# library("slurmjobs")
# job_single(
#   name = "06_qc_scater_scran_metrics", memory = "30G", cores = 1, create_shell = TRUE,
#   task_num = 10
# )

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
