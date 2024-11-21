########################################################################
## Calculate standard quality control metrics based on the singleCellExperiment assay in Seurat's GEX assays
## CSC. Nov-2024
########################################################################

library("SingleCellExperiment")
library("Seurat")
library("scuttle")
library("here")
library("ggplot2")
library("scater")
library("scran")
library("scry")
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
plotDir_reanalyze <- here("plots", "01_preprocessing_QC", "cellrangerARC_reanalyze", "plots_by_sample")

## Load raw data
#unfiltered_path <- here(cellrangerDir_reanalyze, "raw_feature_bc_matrix.h5")

# Check processed_data and plot directories exists
if (!dir.exists(csvDir_reanalyze)) { dir.create(csvDir_reanalyze) }
if (!dir.exists(plotDir_reanalyze)) { dir.create(plotDir_reanalyze) }

## Read the raw_feature_bc_matrix.h5
message("Reading raw feature barcode data corresponding to sample: ", Seurat_base_name) # ../cellrangerARC/S1_Hb_KDM/outs/raw_feature_bc_matrix.h5

set.seed(11112024)


################ Calculate Outliers on GEX multiome assays (By sample) ##############################################

## Load sce.out from cellrangerARC_reanalyze, then transform to sce object
SeuratOBJ <- readRDS(here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", paste0(Seurat_base_name, "_reanalysis.rds")))
# head(sce@meta.data["orig.ident"])
sce <- as.SingleCellExperiment(SeuratOBJ, assay = "RNA")

total_unfiltered_cells <- ncol(sce) # cells in cols

####### Quality control, check low quality cells #######
  
## High mito
is.mito <- grep("MT-", rownames(sce))
sce <- scuttle::addPerCellQC(
  sce,
  subsets = list(Mito = is.mito),
  BPPARAM = BiocParallel::MulticoreParam(4)
)
#colnames(colData(sce))

sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher") # batch = sce$Sample, running one sample at time
## cells pass
table(sce$high_mito)

## low library size
sce$low_sum <- isOutlier(sce$sum, log = TRUE, type = "lower") # , batch = sce$Sample
table(sce$low_sum)

## low detected features
sce$low_detected <- isOutlier(sce$detected, log = TRUE, type = "lower") #, batch = sce$Sample
table(sce$low_detected)

## All low sum are also low detected
table(sce$low_sum, sce$low_detected)

## Annotate cells to remove
sce$discard_auto <- sce$high_mito | sce$low_sum | sce$low_detected
table(sce$discard_auto)

# Filter cells that PASS Outliers and save barcodes filtered
sce_bc_gex <- colnames(sce)[!sce$discard_auto]
total_filtered_cells <-  length(sce_bc)
  
message("Total cells GEX filtered (PASS) from sample ", Seurat_base_name, ": ", total_filtered_cells, " from ", total_unfiltered_cells)

#colnames(sce)[sce$high_mito]
csv_name <- here(csvDir_reanalyze, paste0(Seurat_base_name, "_bc_PASS_isOutliers.csv"))
# NOTE. Uncommend line below if you wish to re-run de outliers barcode-detection and replace the previous csv barcode files
#write.csv(sce_bc_gex, csv_name, row.names=FALSE)

message("Saved ", total_filtered_cells," valid (PASS) barcodes for sample ", Seurat_base_name)

## Build title labels
out_detected <- total_unfiltered_cells - total_filtered_cells
out_detected_p <- round( ((out_detected*100) / total_unfiltered_cells), digits = 2 ) 
caption_label <- paste0(out_detected, " cells (", out_detected_p, "%) with outliers detected. Filtered ", total_filtered_cells, " from ", total_unfiltered_cells)

## Build plot with mito, umi and feature outliers
plot_grid_GEX <- gridExtra::grid.arrange(
  plotColData(sce, x = "orig.ident", y = "subsets_Mito_percent", colour_by = "high_mito", point_size = 1.5) + ggtitle("Mitochondrial percentage") +
    xlab("Sample ID") + ylab("Mito percent"), 
  ## low sum
  plotColData(sce, x = "orig.ident", y = "sum", colour_by = "low_sum", point_size = 1.5) + scale_y_log10() + ggtitle("Total count") +
    xlab("Sample ID") + ylab("Sum UMIs"),    
  ## low genes
  plotColData(sce, x = "orig.ident", y = "detected", colour_by = "low_detected", point_size = 1.5) + scale_y_log10() + ggtitle("Total genes") +
    xlab("Sample ID") + ylab("Sum genes"), 
  ncol = 3,
  top = paste0("Outliers detected. Sample ", Seurat_base_name),
  bottom = caption_label
)

# Save the plot
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_isOutliers_per_sample.png"))
ggsave(filename = plotName, plot = plot_grid_GEX)

plt_GEX_low_sum <- plotColData(sce, x = "orig.ident", y = "sum", colour_by = "low_sum", point_size = 1.5) + scale_y_log10() + ggtitle("Total count") +
  xlab("Sample ID") + ylab("Sum UMIs")

message("Saved plot from isOutliers GEX assay for sample ", Seurat_base_name)
  
# # Mito rate vs n detected features
# plotColData(sce,
#             x = "detected", y = "subsets_Mito_percent",
#             colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5)
# # Detected features vs total count
# plotColData(sce,
#             x = "sum", y = "detected",
#             colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5)
  




################ Calculate Outliers on ATAC multiome assays (By sample) ##############################################

## Load sce.out from cellrangerARC_reanalyze, then transform to sce object

sce <- as.SingleCellExperiment(SeuratOBJ, assay = "ATAC")

total_unfiltered_cells <- ncol(sce) # cells in cols

####### Quality control, check low quality cells #######

## High mito
#is.mito <- grep("MT-", rownames(sce))
sce <- scuttle::addPerCellQC(
  sce,
  #subsets = list(Mito = is.mito),
  BPPARAM = BiocParallel::MulticoreParam(4)
)
#colnames(colData(sce))

table(sce$nCount_ATAC)

## low counts ATAC
sce$low_count_ATAC <- isOutlier(sce$nCount_ATAC, log = TRUE, nmads = 3, type = "lower") # batch = sce$Sample, running one sample at time
table(sce$low_count_ATAC)

## low fragments ATAC
sce$low_feature_ATAC <- isOutlier(sce$nFeature_ATAC, log = TRUE, nmads = 3, type = "lower") # batch = sce$Sample, running one sample at time
table(sce$low_feature_ATAC)

## high NS ATAC
sce$high_NS <- isOutlier(sce$nucleosome_signal, nmads = 3, type = "higher") # batch = sce$Sample, running one sample at time
table(sce$low_count_ATAC)

# # All low sum are also low detected
# table(sce$low_sum, sce$low_detected)

## Annotate cells to remove
sce$discard_auto_atac <- sce$low_count_ATAC | sce$low_feature_ATAC | sce$high_NS
table(sce$discard_auto_atac)

# Filter cells that PASS Outliers and save barcodes filtered
sce_bc_atac <- colnames(sce)[!sce$discard_auto_atac]
total_filtered_cells <-  length(sce_bc)

message("Total cells ATAC filtered (PASS) from sample ", Seurat_base_name, ": ", total_filtered_cells, " from ", total_unfiltered_cells)

#colnames(sce)[sce$high_mito]
csv_name <- here(csvDir_reanalyze, paste0(Seurat_base_name, "_bc_PASS_ATAC_isOutliers.csv"))
#write.csv(sce_bc$Barcode, csv_name, row.names=FALSE)
write.csv(sce_bc_atac, csv_name, row.names=FALSE)

message("Saved ", total_filtered_cells," valid (PASS) barcodes for sample ", Seurat_base_name)

## Build title labels
out_detected <- total_unfiltered_cells - total_filtered_cells
out_detected_p <- round( ((out_detected*100) / total_unfiltered_cells), digits = 2 ) 
caption_label <- paste0(out_detected, " cells (", out_detected_p, "%) with outliers detected. Filtered ", total_filtered_cells, " from ", total_unfiltered_cells)

plt_ATAC_low_sum <- plotColData(sce, x = "orig.ident", y = "nCount_ATAC", colour_by = "low_count_ATAC", point_size = 1.5) + 
  scale_y_log10() + ggtitle("Low-count ATAC") + xlab("Sample ID") + ylab("nCount_ATAC")

# ## Build plot with mito, umi and feature outliers
# plot_grid <- gridExtra::grid.arrange(
#   plotColData(sce, x = "orig.ident", y = "nCount_ATAC", colour_by = "low_count_ATAC") + scale_y_log10() + ggtitle("Low-count ATAC") +
#     xlab("Sample ID") + ylab("nCount_ATAC"), 
#   ## low sum
#   plotColData(sce, x = "orig.ident", y = "nFeature_ATAC", colour_by = "low_feature_ATAC") + scale_y_log10() + ggtitle("Low-feature ATAC") +
#     xlab("Sample ID") + ylab("nFeature_ATAC"),    
#   ## low genes
#   plotColData(sce, x = "orig.ident", y = "nucleosome_signal", colour_by = "high_NS") + scale_y_log10() + ggtitle("High NS") +
#     xlab("Sample ID") + ylab("high_NS"), 
#   ncol = 3,
#   top = paste0("Outliers detected. Sample ", Seurat_base_name),
#   bottom = caption_label
# )


message("Saved plot from isOutliers ATAC assay for sample ", Seurat_base_name)





message("Done!")


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
