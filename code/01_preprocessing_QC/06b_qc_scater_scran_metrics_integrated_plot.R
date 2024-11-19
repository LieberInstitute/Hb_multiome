########################################################################
## Plot quality control metrics based on the multiome `Cell RangerARC-reanalyze`  GEX assays
## CSC. Nov-2024
##
########################################################################

library("SingleCellExperiment")
library("Seurat")
library("scuttle")
library("here")
library("ggplot2")
library("scater")
library("scran")
library("scry")
#library("DropletUtils")
library("gridExtra")
library("stringr")
#library("EnsDb.Hsapiens.v86")
library("sessioninfo")


here::here()

## Read directories

#Seurat_base_name <- "seurat.norm_counts_Harmony_All_subset_Outliers"
#cellrangerDir_reanalyze <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze_outliers", paste0(Seurat_base_name, ".rds"))
Seurat_base_name <- "seurat.norm_counts_Harmony_All"
cellrangerDir_reanalyze <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze", paste0(Seurat_base_name, ".rds"))
plotDir_reanalyze <- here("plots", "01_preprocessing_QC", "cellrangerARC_reanalyze")
#processed-data/03_pseudobulking/cellrangerARC_reanalyze/seurat.norm_counts_Harmony_All.rds

# Check processed_data and plot directories exists
if (!dir.exists(plotDir_reanalyze)) { dir.create(plotDir_reanalyze) }

set.seed(19112024)

## Load Seurat integrated `cellrangerARC_reanalyze` dataset, then transform to sce object

sce.out2 <- readRDS(cellrangerDir_reanalyze)

## Verify data
print(table(sce.out2$orig.ident))
sum(table(sce.out2$orig.ident))

## transform to sce
sce.out2 <- as.SingleCellExperiment(sce.out2, assay = "RNA")
# colnames(colData(sce.out2))

# ## Select those cells detected from custom Cell RangerARC reanalyze pipeline
# v_reanalyze_cells <- Cells(sce.out2)
# 
# ## Subset to keep only valid barcodes from Cell RangerARC reanalyze
# sce <- sce[,sce$Barcode %in% c(v_reanalyze_cells)]

total_unfiltered_cells <- ncol(sce.out2) # cells filtered cells

message("Total cells in the dataset: ", total_unfiltered_cells)

sce <- sce.out2

####### Quality control, check low quality cells #######

## High mito
is.mito <- grep("MT-", rownames(sce))
sce <- scuttle::addPerCellQC(
  sce,
  subsets = list(Mito = is.mito),
  BPPARAM = BiocParallel::MulticoreParam(4)
)
# colnames(colData(sce))

## High mito
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$orig.ident) 
table(sce$high_mito)

## low library size
sce$low_sum <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$orig.ident)
table(sce$low_sum)

## low detected features
sce$low_detected <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$orig.ident)
table(sce$low_detected)

## All low sum are also low detected
table(sce$high_mito, sce$low_detected)
table(sce$low_sum, sce$low_detected)

## Annotate cells to remove
sce$discard_auto <- sce$high_mito | sce$low_sum | sce$low_detected

table(sce$discard_auto)
# FALSE  TRUE
# 4523   527

## discard 9% of nuc
100 * sum(sce$discard_auto) / ncol(sce)
# [1] 9.599

(qc_t <- addmargins(table(sce$orig.ident, sce$discard_auto)))

# Filter/subset cells that PASS Outliers and save barcodes filtered
# sce_bc <- sce[,!sce$discard_auto]
# ncol(sce_bc)
# total_filtered_cells <-  length(sce_bc$discard_auto[sce$discard_auto==FALSE])

total_filtered_cells <- length(sce$discard_auto[sce$discard_auto==FALSE])

message("Total cells with no outliers on `Cell RangerARC-reanalyze` dataset: ", total_filtered_cells, " from ", total_unfiltered_cells)

## Build title labels
out_detected <- total_unfiltered_cells - total_filtered_cells
out_detected_p <- round( ((out_detected*100) / total_unfiltered_cells), digits = 2 )
caption_label <- paste0(out_detected, " cells (", out_detected_p, "%) with outliers detected. Filtered ", total_filtered_cells, " from ", total_unfiltered_cells)

# colnames(colData(sce))

## Save integrated plot by specifi feature 
caption_label <- paste("*Cells to discard: ", length(sce$discard_auto[sce$high_mito==TRUE]), " from ", total_unfiltered_cells)
plt1 <- plotColData(sce, x = "orig.ident", y = "subsets_Mito_percent", colour_by = "high_mito") + ggtitle("Mitochondrial percentage") +
  scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)) +
  xlab("Sample ID") + ylab("Mito percent") +
  labs(caption = caption_label) 
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_high_mito.png"))
ggsave(filename = plotName, plot = plt1, width = 10, height = 5, bg="white")

caption_label <- paste("*Cells to discard: ", length(sce$discard_auto[sce$low_sum==TRUE]), " from ", total_unfiltered_cells)
plt1 <- plotColData(sce, x = "orig.ident", y = "sum", colour_by = "low_sum") + scale_y_log10() + ggtitle("Total count")  +
  scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)) +
  xlab("Sample ID") + ylab("Sum UMIs") +
  labs(caption = caption_label) 
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_low_sum.png"))
ggsave(filename = plotName, plot = plt1, width = 10, height = 5, bg="white")

plot_grid <- gridExtra::grid.arrange(  
  plotColData(sce, x = "orig.ident", y = "subsets_Mito_percent", colour_by = "high_mito") + ggtitle("Mitochondrial percentage") +
    xlab("Sample ID") + ylab("Mito percent") +
    theme(axis.text.x=element_blank()),
  ## low sum
  plotColData(sce, x = "orig.ident", y = "sum", colour_by = "low_sum") + scale_y_log10() + ggtitle("Total count") +
    xlab("Sample ID") + ylab("Sum UMIs") +
    scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)), 
  # ## low genes
  # plotColData(sce, x = "orig.ident", y = "detected", colour_by = "low_detected") + scale_y_log10() + ggtitle("Detected features") +
  #   scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)), 
  nrow = 2,
  top = paste0("Outliers detected. Sample on Cell RangerARC reanalyze"),
  bottom = caption_label
)

# Save the plot
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_isOutliers_ALL.png"))
ggsave(filename = plotName, plot = plot_grid, width = 10, height = 8, bg="white")

print("Saved plots with outliers!")

message("Done!")


# Add job array
library("slurmjobs")
job_single(
  name = "06b_qc_scater_scran_metrics_integrated_plot", memory = "30G", cores = 1, create_shell = TRUE,
  task_num = 10
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
