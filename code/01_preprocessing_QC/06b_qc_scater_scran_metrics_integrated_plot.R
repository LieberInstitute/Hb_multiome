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
library("VennDiagram")
library("gridExtra")
library("stringr")
library("sessioninfo")


here::here()

## Read directories
Seurat_base_name <- "seurat.norm_counts_Harmony_All"
cellrangerDir_reanalyze <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze", paste0(Seurat_base_name, ".rds"))
plotDir_reanalyze <- here("plots", "01_preprocessing_QC", "cellrangerARC_reanalyze")

# Check processed_data and plot directories exists
if (!dir.exists(plotDir_reanalyze)) { dir.create(plotDir_reanalyze) }

set.seed(19112024)

## Load Seurat integrated `Cell RangerARC-reanalyze` dataset, then transform to sce object to calculate outliers

sce <- readRDS(cellrangerDir_reanalyze)

## Verify data
print(table(sce$orig.ident))
sum(table(sce$orig.ident))

## transform to sce
sce <- as.SingleCellExperiment(sce, assay = "RNA")
# colnames(colData(sce))

total_unfiltered_cells <- ncol(sce) # cells filtered cells

message("Total cells in the dataset: ", total_unfiltered_cells)


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

## discard 9% of nuc
total_p <- 100 * sum(sce$discard_auto) / ncol(sce)
# [1] 9.599

qc_t <- addmargins(table(sce$orig.ident, sce$discard_auto))
qc_p <- round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)

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

## Venn diagrams to cross the 3 datasets
set_high_mito <- colnames(sce)[sce$high_mito]
l_hm <- length(set_high_mito)
set_low_umi <- colnames(sce)[sce$low_sum]
l_lsum <- length(set_low_umi)
set_low_detected <- colnames(sce)[sce$low_detected]
l_lgene <- length(set_low_detected)

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_GEX_venn_diagram.png"))
venn.diagram(
  x = list(v_high_mito, v_low_umi, v_low_detected),
  category.names = c(paste0("high_mito (", l_hm, ")"), 
                     paste0("low_umi (",l_lsum, ")"), 
                     paste0("low_genes (", l_lgene, ")")),
  filename = plotName,
  output = FALSE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#0000FF'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#0000FF'),
  rotation = 1,
  main = "Outliers detected",
  sub = "GEX - Cell RangerARC-reanalyze",
  main.cex = 0.7,
  sub.cex = 0.4
)


## Save violin plots integrated plots per metric
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

caption_label <- paste("*Cells to discard: ", length(sce$discard_auto[sce$low_detected==TRUE]), " from ", total_unfiltered_cells)
plt1 <- plotColData(sce, x = "orig.ident", y = "sum", colour_by = "low_detected") + scale_y_log10() + ggtitle("Total genes")  +
  scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)) +
  xlab("Sample ID") + ylab("Sum genes") +
  labs(caption = caption_label) 
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_low_gene.png"))
ggsave(filename = plotName, plot = plt1, width = 10, height = 5, bg="white")

## saved plot with all metrics 
caption_label <- paste0("*Cells to discard: ", length(sce$discard_auto[sce$discard_auto==TRUE]), " (", round(total_p, digits = 2) ,"%) from ", total_unfiltered_cells)
plot_grid <- gridExtra::grid.arrange(  
  plotColData(sce, x = "orig.ident", y = "subsets_Mito_percent", colour_by = "high_mito") + ggtitle("Mitochondrial percentage") +
    xlab("Sample ID") + ylab("Mito percent") +
    theme(axis.text.x=element_blank()),
  ## low sum
  plotColData(sce, x = "orig.ident", y = "sum", colour_by = "low_sum") + scale_y_log10() + ggtitle("Total count") +
    xlab("Sample ID") + ylab("Sum UMIs") +
    theme(axis.text.x=element_blank()),
  ## low genes
  plotColData(sce, x = "orig.ident", y = "sum", colour_by = "low_detected") + scale_y_log10() + ggtitle("Total feature")  +
    scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)) +
    xlab("Sample ID") + ylab("Sum genes") +
    scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)), 
  nrow = 3,
  top = paste0("Outliers detected on `Cell RangerARC-reanalyze` dataset"),
  bottom = caption_label
)

# Save the plot
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_all_isOutliers.png"))
ggsave(filename = plotName, plot = plot_grid, width = 10, height = 12, bg="white")

print("Saved plots with outliers!")

message("Done!")


# # Add job array
# library("slurmjobs")
# job_single(
#   name = "06b_qc_scater_scran_metrics_integrated_plot", memory = "30G", cores = 1, create_shell = TRUE,
#   task_num = 10
# )

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
