########################################################################
## Plot Coverage Plots on WNN clusters
##
## Authors. CSC
## Date. April 04, 2025
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## Note. Seurat objects were created with module load conda_R/4.3.x
########################################################################

library("Seurat")
library("Signac")
library("ggplot2")
library("tidyverse")
library("here")

## input directories

here()

# Check/create directories

## clusters renamed for Spatial-Registration on Visium project
inputRDS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "08_wnn_gene_expression_plts_renamed_idents"
)
inputCVS_Dir <- here(
  "processed-data",
  "05_Clustering_ARCr",
  "02_Hb_celltypes_from_seurat_reanalyze_v3",
  "cvs_files_markers"
)
plotDir <- here(
  "plots",
  "06_peak_calling"
)

## Check directories
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}

## Load Seurat
# Use Seurat with clusters renamed for Spatial-Registration on Visium project
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
levels(SeuratOBJ)

DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])
Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+")) 
# C.leiden_lsi_r2_renamed_visium


message("Link Peaks to Genes (cis-regulatory analysis)")

## Find peaks that are correlated with the expression of nearby genes
## - link peaks to gene expression using correlation 
## - For each gene, LinkPeaks() computes the correlation coefficient (CC) between the gene expression and accessibility of each peak within a given distance from the gene TSS, and computes an expected CC for each peak given the GC content, accessibility, and length of the peak. The expected coefficient values for the peak are then used to compute a z-score and p-value.

atac <- LinkPeaks(
  object = SeuratOBJ,
  peak.assay = "peaks",
  expression.assay = "RNA",  # Make sure RNA assay is integrated
  genes.use = NULL,          # Or supply vector of gene names if you're interested in a subset
  method = "pearson",         # I am starting with default settings
  distance = 1e5             # Cis distance (e.g., 100kb window)
)

## Optionally accounting for covariates

# Coverage and gene link plot for a gene of interest
CoveragePlot(
  object = atac,
  region = "GENE_NAME",      # Replace with e.g. "CD14"
  features = "GENE_NAME",
  expression.assay = "RNA",
  extend.upstream = 10000,
  extend.downstream = 10000
)
