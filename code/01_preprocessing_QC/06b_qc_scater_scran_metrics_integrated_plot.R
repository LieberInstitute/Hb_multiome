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
library("RColorBrewer")
library("gridExtra")
library("stringr")
library("sessioninfo")


here::here()
# Prepare a palette of 3 colors with R colorbrewer:
myCol <- brewer.pal(3, "Pastel2")


## Read directories
Seurat_base_name <- "seurat.norm_counts_Harmony_All"
cellrangerDir_reanalyze <- here("processed-data", "03_pseudobulking", "cellrangerARC_reanalyze", paste0(Seurat_base_name, ".rds"))
plotDir_reanalyze <- here("plots", "01_preprocessing_QC", "cellrangerARC_reanalyze")

# Check processed_data and plot directories exists
if (!dir.exists(plotDir_reanalyze)) { dir.create(plotDir_reanalyze) }

set.seed(19112024)

## Load Seurat integrated `Cell RangerARC-reanalyze` dataset, then transform to sce object to calculate outliers

SeuratOBJ <- readRDS(cellrangerDir_reanalyze)

## Verify data
print(table(SeuratOBJ$orig.ident))
sum(table(SeuratOBJ$orig.ident))
# [1] 62950

################ Prepare plots to visualize outliers on GEX multiome assays ###################

## transform to sce
sce <- as.SingleCellExperiment(SeuratOBJ, assay = "RNA")
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

## low detected features/genes
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

message("Cells without outliers: ", total_filtered_cells, " from ", total_unfiltered_cells)

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


## Prepare and plot Venn diagram with GEX high mito, low-umi and low-feature detected
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_GEX_VENNd_outliers_detected.png"))
venn.diagram(
  x = list(set_high_mito, set_low_umi, set_low_detected),
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

# ## Detected mito vs total count
# # Mito rate vs n detected features
# plotColData(sce,
#             x = "detected", y = "subsets_Mito_percent",
#             colour_by = "discard_auto", point_size = 2, point_alpha = 0.5
# )
# # Detected features vs total count
# plotColData(sce,
#             x = "sum", y = "detected",
#             colour_by = "discard_auto", point_size = 2, point_alpha = 0.5
# )


## Save violin plots with all samples integrated per metric
caption_label <- paste("*Cells to discard: ", length(sce$discard_auto[sce$high_mito]), " from ", total_unfiltered_cells)
plt_hm <- plotColData(sce, x = "orig.ident", y = "subsets_Mito_percent", colour_by = "high_mito", point_size = 0.5) + 
  ggtitle("Mitochondrial percentage") +
  scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)) +
  xlab("Sample ID") + ylab("Mito percent") +
  labs(caption = caption_label) 
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_GEX_OUTLIERS_high_mito.png"))
ggsave(filename = plotName, plot = plt_hm, width = 10, height = 5, bg="white")

caption_label <- paste("*Cells to discard: ", length(sce$discard_auto[sce$low_sum]), " from ", total_unfiltered_cells)
plt_low_sum_gex <- plotColData(sce, x = "orig.ident", y = "sum", colour_by = "low_sum", point_size = 0.5) + scale_y_log10() +
  ggtitle("Total count")  +
  scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)) +
  xlab("Sample ID") + ylab("Sum UMIs") +
  labs(caption = caption_label) 
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_GEX_OUTLIERS_low_sum.png"))
ggsave(filename = plotName, plot = plt_low_sum_gex, width = 10, height = 5, bg="white")

caption_label <- paste("*Cells to discard: ", length(sce$discard_auto[sce$low_detected]), " from ", total_unfiltered_cells)
plt_low_genes_gex <- plotColData(sce, x = "orig.ident", y = "sum", colour_by = "low_detected", point_size = 0.5) + scale_y_log10() + 
  ggtitle("Total genes")  +
  scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)) +
  xlab("Sample ID") + ylab("Sum genes") +
  labs(caption = caption_label) 
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_GEX_OUTLIERS_low_gene.png"))
ggsave(filename = plotName, plot = plt_low_genes_gex, width = 10, height = 5, bg="white")

## saved plot with all metrics 
caption_label <- paste0("*Cells to discard: ", length(sce$discard_auto[sce$discard_auto]), " (", round(total_p, digits = 2) ,"%) from ", total_unfiltered_cells)
plot_grid <- gridExtra::grid.arrange(  
  plt_hm + ggtitle("Mitochondrial percentage") + labs(caption = "") +
    theme(axis.text.x=element_blank()),
  plt_low_sum_gex + ggtitle("Total count") + labs(caption = "") +
    theme(axis.text.x=element_blank()),
  plt_low_genes_gex + ggtitle("Total feature") + labs(caption = "") +
    scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)) +
    scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)), 
  nrow = 3,
  top = paste0("Outliers detected on GEX `Cell RangerARC-reanalyze` dataset"),
  bottom = caption_label
)

# Save the plot
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_GEX_ALL_OUTLIERS.png"))
ggsave(filename = plotName, plot = plot_grid, width = 10, height = 12, bg="white")

print("Saved GEX plots with outliers!")

message("GEX Done!")




################ Calculate Outliers on ATAC multiome assays (By sample) ##############################################

total_filtered_cells <- 0
total_unfiltered_cells <- 0

sce_atac <- as.SingleCellExperiment(SeuratOBJ, assay = "ATAC")

total_unfiltered_cells <- ncol(sce_atac) # cells in cols

####### Quality control, check low quality cells #######

## Perform quality metrics
sce_atac <- scuttle::addPerCellQC(
  sce_atac,
  BPPARAM = BiocParallel::MulticoreParam(4)
)
#colnames(colData(sce_atac))

#table(sce_atac$nCount_ATAC)

## low counts ATAC
sce_atac$low_sum_ATAC <- isOutlier(sce_atac$nCount_ATAC, log = TRUE, nmads = 3, type = "lower", batch = sce$orig.ident) 
table(sce_atac$low_sum_ATAC)

## low fragments ATAC
sce_atac$low_feature_ATAC <- isOutlier(sce_atac$nFeature_ATAC, log = TRUE, nmads = 3, type = "lower", batch = sce$orig.ident)
#table(sce_atac$low_feature_ATAC)

## high NS ATAC
max(sce_atac$nucleosome_signal); min(sce_atac$nucleosome_signal)
sce_atac$high_NS <- isOutlier(sce_atac$nucleosome_signal, nmads = 3, type = "higher", batch = sce$orig.ident)
table(sce_atac$high_NS)
#sce_atac$nucleosome_signal[sce_atac$high_NS]

## Low TSS.enrichment
max(sce_atac$TSS.enrichment); min(sce_atac$TSS.enrichment)
sce_atac$low_TSS <- isOutlier(sce_atac$TSS.enrichment, nmads = 3, type = "lower") # batch = sce_atac$Sample, running one sample at time
table(sce_atac$low_TSS)
sce$TSS.enrichment[sce_atac$low_TSS]
#[1] 0.9571824 1.1462222 1.1901142 1.0215137 1.1916909
#pmatch(colnames(sce_atac)[sce_atac$high_NS], colnames(sce_atac)[sce_atac$low_TSS])
## Annotate cells to remove
sce_atac$discard_auto_atac <- sce_atac$low_sum_ATAC | sce_atac$low_feature_ATAC | sce_atac$low_TSS #sce_atac$high_NS
table(sce_atac$discard_auto_atac)

## discard 9% of nuc
total_p <- 100 * sum(sce_atac$discard_auto_atac) / ncol(sce_atac)
#[1] 2.27641

qc_t <- addmargins(table(sce_atac$orig.ident, sce_atac$discard_auto_atac))
qc_p <- round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)

# Filter/subset cells that PASS Outliers and save barcodes filtered
# sce_bc <- sce[,!sce$discard_auto]
# ncol(sce_bc)
# total_filtered_cells <-  length(sce_bc$discard_auto[sce$discard_auto==FALSE])

total_filtered_cells <- length(sce_atac$discard_auto_atac[sce_atac$discard_auto_atac==FALSE])

message("Cells without outliers: ", total_filtered_cells, " from ", total_unfiltered_cells)

## Build title labels
out_detected <- total_unfiltered_cells - total_filtered_cells
out_detected_p <- round( ((out_detected*100) / total_unfiltered_cells), digits = 2 )
caption_label <- paste0(out_detected, " cells (", out_detected_p, "%) with outliers detected. Filtered ", total_filtered_cells, " from ", total_unfiltered_cells)

# colnames(colData(sce))

## Venn diagrams to cross the 3 datasets
set_low_atac <- colnames(sce_atac)[sce_atac$low_sum_ATAC]
l_latac <- length(set_low_atac)
set_low_detected_atac <- colnames(sce_atac)[sce_atac$low_feature_ATAC]
l_lsum_atac <- length(set_low_detected_atac)
set_low_TSS <- colnames(sce_atac)[sce_atac$low_TSS]
l_lTSS <- length(set_low_TSS)


## Prepare and plot Venn diagram with GEX high mito, low-umi and low-feature detected
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_ATAC_VENNd_outliers_detected.png"))
venn.diagram(
  x = list(set_low_atac, set_low_detected_atac, set_low_TSS),
  category.names = c(paste0("nCount_ATAC (", l_latac, ")"), 
                     paste0("nFeature_ATAC (",l_lsum_atac, ")"), 
                     paste0("low_TSS (", l_lTSS, ")")),
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
  sub = "ATAC - Cell RangerARC-reanalyze",
  main.cex = 0.7,
  sub.cex = 0.4
)

#rm seurat.norm_counts_Harmony_All_ATAC_VENNd_outliers_detected.png.2024-11-21_21-27-16.933048.log

## Save violin plots with all samples integrated per metric
caption_label <- paste("*Cells to discard: ", length(sce$discard_auto[sce$high_mito]), " from ", total_unfiltered_cells)
plt_hm <- plotColData(sce, x = "orig.ident", y = "subsets_Mito_percent", colour_by = "high_mito", point_size = 0.5) + 
  ggtitle("Mitochondrial percentage") +
  scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)) +
  xlab("Sample ID") + ylab("Mito percent") +
  labs(caption = caption_label) 
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_GEX_OUTLIERS_high_mito.png"))
ggsave(filename = plotName, plot = plt_hm, width = 10, height = 5, bg="white")

caption_label <- paste("*Cells to discard: ", length(sce$discard_auto[sce$low_sum]), " from ", total_unfiltered_cells)
plt_low_sum_gex <- plotColData(sce, x = "orig.ident", y = "sum", colour_by = "low_sum", point_size = 0.5) + scale_y_log10() +
  ggtitle("Total count")  +
  scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)) +
  xlab("Sample ID") + ylab("Sum UMIs") +
  labs(caption = caption_label) 
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_GEX_OUTLIERS_low_sum.png"))
ggsave(filename = plotName, plot = plt_low_sum_gex, width = 10, height = 5, bg="white")

caption_label <- paste("*Cells to discard: ", length(sce$discard_auto[sce$low_detected]), " from ", total_unfiltered_cells)
plt_low_genes_gex <- plotColData(sce, x = "orig.ident", y = "sum", colour_by = "low_detected", point_size = 0.5) + scale_y_log10() + 
  ggtitle("Total genes")  +
  scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)) +
  xlab("Sample ID") + ylab("Sum genes") +
  labs(caption = caption_label) 
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_GEX_OUTLIERS_low_gene.png"))
ggsave(filename = plotName, plot = plt_low_genes_gex, width = 10, height = 5, bg="white")

## saved plot with all metrics 
caption_label <- paste0("*Cells to discard: ", length(sce$discard_auto[sce$discard_auto]), " (", round(total_p, digits = 2) ,"%) from ", total_unfiltered_cells)
plot_grid <- gridExtra::grid.arrange(  
  plt_hm + ggtitle("Mitochondrial percentage") + labs(caption = "") +
    theme(axis.text.x=element_blank()),
  plt_low_sum_gex + ggtitle("Total count") + labs(caption = "") +
    theme(axis.text.x=element_blank()),
  plt_low_genes_gex + ggtitle("Total feature") + labs(caption = "") +
    scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)) +
    scale_x_discrete(labels = ~ str_wrap(gsub('_', ' ', .x), 10)), 
  nrow = 3,
  top = paste0("Outliers detected on GEX `Cell RangerARC-reanalyze` dataset"),
  bottom = caption_label
)

# Save the plot
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_GEX_ALL_OUTLIERS.png"))
ggsave(filename = plotName, plot = plot_grid, width = 10, height = 12, bg="white")

print("Saved GEX plots with outliers!")

message("ATAC Done!")


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
