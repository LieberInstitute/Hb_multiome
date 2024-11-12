########################################################################
## Calculate standard quality control metrics based on the singleCellExperiment assay in Seurat's ATAC assays
## CSC. Nov-2024
########################################################################

library("SingleCellExperiment")
library("Seurat")
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

## Scans arguments invoked from slurm job shell sh
sample_tmp <- commandArgs(trailingOnly = TRUE)
# For testing:
# sample_tmp <- "S3_Hb_KDM_reanalysis, S3_Hb_KDM"
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

set.seed(12112024)


# ############### FINDING HIGH MITO ##############################################

## Call function
process_sample <- function(sample_path, sce_out_path, fdr_threshold = 0.001) {
  # fdr_threshold = 0.001   # Only need it if you are processing empyDrops results
  # sample_path <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/cellrangerARC/S3_Hb_KDM/outs/raw_feature_bc_matrix.h5"
  
  # Load the sample data
  sce <- read10xCounts(sample_path)
  # class: SingleCellExperiment 
  # dim: 181201 702703 
  
  # Unifying feature names
  rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
  
  # Map the IDs
  location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, column="SEQNAME", keytype="GENEID")
  
  ## Load sce.out from cellrangerARC_reanalyze Dir
  SeuratOBJ <- readRDS(here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", paste0(Seurat_base_name, "_reanalysis.rds")))
  DefaultAssay(SeuratOBJ) <- "ATAC"
  
  #table(SeuratOBJ$nucleosome_signal)
  ## Add new meta data with NS position
  SeuratOBJ$nucleosome_signal_p <- round(SeuratOBJ$nucleosome_signal, digits = 0)
  table(SeuratOBJ$nucleosome_signal_p) 
  sum(SeuratOBJ$nucleosome_signal_p == 0)
  ## Add new meta-data with NS positions
  SeuratOBJ$nucleosome_position_pos <- ifelse(SeuratOBJ$nucleosome_signal_p == 0, "NucleosomeFree", 
                                              ifelse(SeuratOBJ$nucleosome_signal_p == 2, "Dinucleosomes", "Multinucleosomes")) 
  #table(SeuratOBJ$nucleosome_position_pos)
  #sum(SeuratOBJ$nucleosome_position_pos == "NucleosomeFree")
  sum(!SeuratOBJ$nucleosome_position_pos == "NucleosomeFree")
  
  
  
  ## Convert ATAC assay from Seurat to singleCellexperiment
  
  sce.out2 <- as.SingleCellExperiment(SeuratObj, assay = "ATAC")
  ## Select those cells detected from custom Cell RangerARC reanalyze pipeline
  v_reanalyze_cells <- Cells(sce.out2)
  #dim(sce) ## [1] 181201 702703
  #colData(sce)
  sce <- sce[,sce$Barcode %in% c(v_reanalyze_cells)]
  #dim(sce) #[1] 181201   5050
  unfiltered <- sce
  
  total_unfiltered_cells <- dim(assay(unfiltered))[2]
  
  # Quality control
  # Filtering on the mitochondrial proportion
  is.mito <- grep("MT-", rownames(sce))
  stats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
  
  ## For reference, these are result derived from reanalyze (5050 cells) vs emptyDrops for sample S03 (5802 cells)
  #summary(stats$subsets_Mito_percent) 
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.0000  0.4765  1.0858  2.5113  2.1428 84.0290
  
  ## For comparison purposes: result derived from emptyDrops for sample S03
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 0.0000  0.5258  1.1725  2.5728  2.3095 71.2544
  #stats$subsets_Mito_percent[is.na(stats$subsets_Mito_percent)] <- 0

  #stats <- perCellQCMetrics(sce, subsets=list(Mito=which(location=="MT")))
  high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
  table(high.mito[high.mito==FALSE])
  
  sce <- sce[,!high.mito]
  total_filtered_cells <- dim(assay(sce))[2]
  
  ## store this in the colData() of our SingleCellExperiment object for future reference
  colData(unfiltered) <- cbind(colData(unfiltered), stats)
  unfiltered$discard <- high.mito
  # colnames(colData(unfiltered))
  
  csv_name <- here(csvDir_reanalyze, paste0(Seurat_base_name, "_isOutliers_valid_barcodes.csv"))
  write.csv(sce$Barcode, csv_name)
  print(paste0("Saved valid (true) barcodes from isOutliers for sample ", Seurat_base_name))

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
  ggsave(filename = plotName, plot = plot_grid)
  print(paste0("Saved plot from isOutliers for sample ", Seurat_base_name))
  
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
