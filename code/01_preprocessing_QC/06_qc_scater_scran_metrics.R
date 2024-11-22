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
library("VennDiagram")
library("RColorBrewer")
library("EnsDb.Hsapiens.v86")
library("sessioninfo")

## Read directories

here::here()
myCol <- brewer.pal(3, "Pastel2")

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
plotDir_reanalyze <- here("plots", "01_preprocessing_QC", "cellrangerARC_reanalyze", "plots_by_sample")

## Load raw data
#unfiltered_path <- here(cellrangerDir_reanalyze, "raw_feature_bc_matrix.h5")

# Check processed_data and plot directories exists
if (!dir.exists(csvDir_reanalyze)) { dir.create(csvDir_reanalyze) }
if (!dir.exists(plotDir_reanalyze)) { dir.create(plotDir_reanalyze) }

## Read the raw_feature_bc_matrix.h5
message("Reading raw feature barcode data corresponding to sample: ", Seurat_base_name) # ../cellrangerARC/S1_Hb_KDM/outs/raw_feature_bc_matrix.h5

set.seed(11112024)


## Load sce.out from cellrangerARC_reanalyze, then transform to sce object

SeuratOBJ <- readRDS(here("processed-data", "01_preprocessing_QC", "cellrangerARC_reanalyze", paste0(Seurat_base_name, "_reanalysis.rds")))
# head(sce@meta.data["orig.ident"])
sce <- as.SingleCellExperiment(SeuratOBJ, assay = "RNA")

total_unfiltered_cells <- ncol(sce) # cells in cols


################ Calculate Outliers on GEX multiome assays (By sample) ##############################################

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
total_filtered_cells <-  length(sce_bc_gex)
  
message("Total cells GEX filtered (PASS) from sample ", Seurat_base_name, ": ", total_filtered_cells, " from ", total_unfiltered_cells)

#colnames(sce)[sce$high_mito]
csv_name <- here(csvDir_reanalyze, paste0(Seurat_base_name, "_bc_PASS_isOutliers.csv"))
# NOTE. Uncommend line below if you wish to re-run de outliers barcode-detection and replace the previous csv barcode files
#write.csv(sce_bc_gex, csv_name, row.names=FALSE)

message("Saved ", total_filtered_cells," valid (PASS) barcodes for sample ", Seurat_base_name)

## Build title labels
out_detected <- total_unfiltered_cells - total_filtered_cells
out_detected_p <- round( ((out_detected*100) / total_unfiltered_cells), digits = 2 ) 
caption_label_GEX <- paste0("GEX: ", out_detected, " (", out_detected_p, "%)", " outliers detected. ", total_filtered_cells, " good quality cells of ", total_unfiltered_cells)

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
  top = paste0("Outliers detected. Sample ", unlist(strsplit(Seurat_base_name, split = "_"))[1]),
  bottom = caption_label_GEX
)

# Save the plot
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_isOutliers_per_sample.png"))
ggsave(filename = plotName, plot = plot_grid_GEX)

## build plot with low-counts for GEX to compare later
plt_GEX_low_sum <- plotColData(sce, x = "orig.ident", y = "sum", colour_by = "low_sum", point_size = 1.5) + scale_y_log10() + ggtitle("Total count") +
  xlab("Sample ID") + ylab("Sum UMIs")  + theme(legend.position="none")


## Venn diagrams to cross the 3 sets of metrics
set_high_mito <- colnames(sce)[sce$high_mito]
l_high_mito <- length(set_high_mito)

set_low_gex <- colnames(sce)[sce$low_sum]
l_lsum_gex <- length(set_low_gex)

set_low_feature_gex <- colnames(sce)[sce$low_detected]
l_lgenes_gex <- length(set_low_feature_gex)

# Venn diagram for ATAC metrics 
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_GEX_VENNd_low_counts.png"))
venn.diagram(
  x = list(set_high_mito, set_low_gex, set_low_feature_gex),
  category.names = c(paste0("high_mito (", l_high_mito, ")"),
                     paste0("low_sum (",l_lsum_gex, ")"), 
                     paste0("low_genes (", l_lgenes_gex, ")")),
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


message("Saved plots with outliers  on GEX assay sample ", Seurat_base_name)
  
# # Mito rate vs n detected features
# plotColData(sce,
#             x = "detected", y = "subsets_Mito_percent",
#             colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5)
# # Detected features vs total count
# plotColData(sce,
#             x = "sum", y = "detected",
#             colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5)
  




################ Calculate Outliers on ATAC multiome assays (By sample) ##############################################

total_filtered_cells <- 0
total_unfiltered_cells <- 0

sce_atac <- as.SingleCellExperiment(SeuratOBJ, assay = "ATAC")

total_unfiltered_cells <- ncol(sce_atac) # cells in cols

####### Quality control, check low quality cells #######

## Perform quality metrics
sce_atac <- scuttle::addPerCellQC(
  sce_atac,
  #subsets = list(Mito = is.mito),
  BPPARAM = BiocParallel::MulticoreParam(4)
)
#colnames(colData(sce_atac))

#table(sce_atac$nCount_ATAC)

## low counts ATAC
sce_atac$low_sum_ATAC <- isOutlier(sce_atac$nCount_ATAC, log = TRUE, nmads = 3, type = "lower") # batch = sce_atac$Sample, running one sample at time
table(sce_atac$low_sum_ATAC)

## low fragments ATAC
sce_atac$low_feature_ATAC <- isOutlier(sce_atac$nFeature_ATAC, log = TRUE, nmads = 3, type = "lower") # batch = sce_atac$Sample, running one sample at time
#table(sce_atac$low_feature_ATAC)

## high NS ATAC
max(sce_atac$nucleosome_signal); min(sce_atac$nucleosome_signal)
sce_atac$high_NS <- isOutlier(sce_atac$nucleosome_signal, nmads = 3, type = "higher") # batch = sce_atac$Sample, running one sample at time
table(sce_atac$high_NS)
sce_atac$nucleosome_signal[sce_atac$high_NS]
#[1] 1.440824 1.396084 1.466782 1.432133 1.540221

## Low TSS.enrichment
max(sce_atac$TSS.enrichment); min(sce_atac$TSS.enrichment)
sce_atac$low_TSS <- isOutlier(sce_atac$TSS.enrichment, nmads = 3, type = "lower") # batch = sce_atac$Sample, running one sample at time
table(sce_atac$low_TSS)
sce$TSS.enrichment[sce_atac$low_TSS]
#[1] 0.9571824 1.1462222 1.1901142 1.0215137 1.1916909
#pmatch(colnames(sce_atac)[sce_atac$high_NS], colnames(sce_atac)[sce_atac$low_TSS])
## Annotate cells to remove
sce_atac$discard_auto_atac <- sce_atac$low_sum_ATAC | sce_atac$low_feature_ATAC | sce_atac$high_NS | sce_atac$low_TSS
table(sce_atac$discard_auto_atac)

## Filter cells that PASS Outliers and save barcodes filtered
sce_bc_atac <- colnames(sce_atac)[!sce_atac$discard_auto_atac]

total_filtered_cells <-  length(sce_bc_atac)

message("Total cells ATAC filtered (PASS) from sample ", Seurat_base_name, ": ", total_filtered_cells, " from ", total_unfiltered_cells)

#colnames(sce_atac)[sce_atac$high_mito]
csv_name <- here(csvDir_reanalyze, paste0(Seurat_base_name, "_bc_PASS_ATAC_isOutliers.csv"))
write.csv(sce_bc_atac, csv_name, row.names=FALSE)

message("Saved ", total_filtered_cells," valid (PASS) barcodes for sample ", Seurat_base_name)

## Build title labels
out_detected <- total_unfiltered_cells - total_filtered_cells
out_detected_p <- round( ((out_detected*100) / total_unfiltered_cells), digits = 2 ) 
caption_label_ATAC <- paste0("ATAC: ", out_detected, " (", out_detected_p, "%) outliers detected. ", total_filtered_cells, " good quality cells of ",
                             total_unfiltered_cells)

## Plot GEX total-counts and ATAC total-counts
plot_grid_multiome <- gridExtra::grid.arrange(
  plt_GEX_low_sum,
  plotColData(sce_atac, x = "orig.ident", y = "nCount_ATAC", colour_by = "low_sum_ATAC", point_size = 1.5) + 
    scale_y_log10() + ggtitle("Low-count ATAC") + xlab("Sample ID") + ylab("nCount_ATAC") + theme(legend.title=element_blank()), 
  ncol = 2,
  top = paste0("Outliers cells GEX and ATAC. Sample ", unlist(strsplit(Seurat_base_name, split = "_"))[1]),
  bottom = paste("Total low-sum on GEX: ", length(colnames(sce)[sce$low_sum]), "\n", "Total low-sum on ATAC: ", length(colnames(sce_atac)[sce_atac$low_sum_ATAC])) 
)

# Save plot comparing total counts in RNA and ATAC
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_isOutliers_GEX_ATAC_total_counts.png"))
ggsave(filename = plotName, plot = plot_grid_multiome, width = 8, height = 10)

# ## Plot NS and TSS ATAC metrics
# plot_grid_multiome <- gridExtra::grid.arrange(
#   plotColData(sce_atac, x = "orig.ident", y = "nucleosome_signal", colour_by = "high_NS", point_size = 1.5) + 
#     scale_y_log10() + ggtitle("High nucleosome signal") + xlab("Sample ID") + ylab("high_NS")  + theme(legend.position="none"),
#   plotColData(sce_atac, x = "orig.ident", y = "TSS.enrichment", colour_by = "low_TSS", point_size = 1.5) + 
#     scale_y_log10() + ggtitle("Low TSS.enrichment") + xlab("Sample ID") + ylab("low_TSS") + theme(legend.title=element_blank()), 
#   ncol = 2,
#   top = paste0("Outliers cells NS and TSS. Sample ", unlist(strsplit(Seurat_base_name, split = "_"))[1]),
#   bottom = paste("High-NS:", length(colnames(sce_atac)[sce_atac$high_NS]), "Low-TSS:", length(colnames(sce_atac)[sce_atac$low_TSS])) 
# )

## Venn diagrams to cross the 3 sets of metrics
set_low_atac <- colnames(sce_atac)[sce_atac$low_sum_ATAC]
l_lsum_atac <- length(set_low_atac)

set_low_feature_atac <- colnames(sce_atac)[sce_atac$low_feature_ATAC]
l_lgenes_atac <- length(set_low_feature_atac)

set_low_TSS <- colnames(sce_atac)[sce_atac$low_TSS]
l_low_TSS <- length(set_low_TSS)

# Venn diagram for ATAC metrics 
plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_ATAC_VENNd_low_counts.png"))
venn.diagram(
  x = list(set_low_atac, set_low_feature_atac, set_low_TSS),
  category.names = c(paste0("low_sum (",l_lsum_atac, ")"), 
                     paste0("low_genes (", l_lgenes_atac, ")"),
                     paste0("low_tss (", l_low_TSS, ")")),
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

# ## Build plot with low counts, features and high ns
# plot_grid <- gridExtra::grid.arrange(
#   plotColData(sce_atac, x = "orig.ident", y = "nCount_ATAC", colour_by = "low_sum_ATAC") + scale_y_log10() + ggtitle("Low-count ATAC") +
#     xlab("Sample ID") + ylab("nCount_ATAC"),
#   ## low sum
#   plotColData(sce_atac, x = "orig.ident", y = "nFeature_ATAC", colour_by = "low_feature_ATAC") + scale_y_log10() + ggtitle("Low-feature ATAC") +
#     xlab("Sample ID") + ylab("nFeature_ATAC"),
#   ## low genes
#   plotColData(sce_atac, x = "orig.ident", y = "nucleosome_signal", colour_by = "high_NS") + scale_y_log10() + ggtitle("High NS") +
#     xlab("Sample ID") + ylab("high_NS"),
#   ncol = 3,
#   top = paste0("Outliers detected. Sample ", Seurat_base_name),
#   bottom = caption_label_ATAC
# )
# 
# # Save plot comparing total counts in RNA and ATAC
# plotName <- here(plotDir_reanalyze, paste0(Seurat_base_name, "_isOutliers_ATAC.png"))
# ggsave(filename = plotName, plot = plot_grid_multiome, width = 8, height = 10)


message("Saved plots with outliers on ATAC assay sample ", Seurat_base_name)





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
