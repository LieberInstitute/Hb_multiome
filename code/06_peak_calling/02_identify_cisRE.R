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
library("BSgenome.Hsapiens.UCSC.hg38")
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
  "17_wnn_clustering_final_ct"
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
# Use Seurat with WNN clusters shared with Visium projects (ge. to plot Spatial-Registration, etc)
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"
seurat_name <- here(inputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(here(inputRDS_Dir, Seurat_base_name))
levels(SeuratOBJ)

DefaultAssay(SeuratOBJ) <- "ATAC"
class(SeuratOBJ[["ATAC"]])
Seurat_base_name <- str_extract(seurat_name, regex("C\\.\\w+")) 
# C.leiden_lsi_r2_renamed_visium

#colnames(SeuratOBJ@meta.data)
message("Processing ", length(unique(SeuratOBJ$seurat_clusters)), " clusters")


# ## ========================================================================== ##
# 
# ## (2) This identifies cis-regulatory elements by linking chromatin-accessible peaks to gene expression using correlation (and optionally accounting for covariates).
# ## Here I am only considering the 13 HABENULA WNN clusters
# 
# ## ========================================================================== ##
# 
# # Find annotated Habenula clusters and extract their top-5 most expressed genes to link peaks
# 
# # First read all DEG
# 
# DEG_file_name <- "WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_cellTypes_integrated_top50.csv"
# DEG_file_name <- here(inputCVS_Dir, DEG_file_name)
# df_cluster_names <- read.csv(DEG_file_name)
# df_cluster_names <- df_cluster_names |> drop_na(cell_type)
# head(df_cluster_names)
# # p_val avg_log2FC pct.1 pct.2 p_val_adj cluster     gene            cell_type
# # 1     0   4.336764 0.938 0.100         0       1 OTX2-AS1        DD_Inhib.Thal
# # 2     0   3.849508 0.900 0.083         0       1      KIT        DD_Inhib.Thal
# # 3     0   3.682782 0.927 0.122         0       1    MEIS2      LB_Thalamus/MDm
# 
# ## Identified and subset clusters annotated for `habenula`.
# 
# message("Cluster-IDs from `WNN`")
# 
# SeuOBJ_clusters <- Idents(SeuratOBJ)
# hb_clusters <- unlist(levels(SeuOBJ_clusters))
# ## Get top 5. Filter habenula clusters only
# no_hb_clust = list()
# for (idx in seq_along(hb_clusters)) {
#   if (nchar(hb_clusters[idx]) <= 4) {
#     no_hb_clust <- append(no_hb_clust, hb_clusters[idx])
#   }
# }
# no_hb_clust <- c(unlist(no_hb_clust))
# hb_clusters <- hb_clusters[!hb_clusters %in% c(no_hb_clust)]
# as.vector(hb_clusters)
# # [1] "C.05 DD_LHb" "C.07 DD_MHb" "C.10 DD_MHb" "C.11 DD_MHb" "C.14 DD_MHb"
# # [6] "C.16 DD_MHb" "C.18 DD_LHb" "C.23 DD_LHb" "C.24 DD_LHb" "C.30 DD_LHb"
# # [11] "C.33 DD_LHb" "C.36 DD_MHb" "C.40 DD_LHb"
# 
# ##  Use length of clusters to extract clusters IDs
# hb_clusters <- as.integer(substr(hb_clusters, 3, 4))
# # [1]  5  7 10 11 14 16 18 23 24 30 33 36 40
# 
# ## filter the top 5 most expressed genes
# 
# unique(df_cluster_names$cluster)
# top5 <- df_cluster_names |>
#   filter(cluster %in% hb_clusters) |>
#   group_by(cluster) |>
#   top_n(n = 5, wt = avg_log2FC) |>
#   select(cluster, gene)
# head(top5)
# # p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene    cell_type
# # <dbl>      <dbl> <dbl> <dbl>     <dbl>   <int> <chr>   <chr>    
# # 1     0       3.04 0.768 0.145         0       5 RFTN1   DD_LHb   
# # 2     0       3.03 0.798 0.194         0       5 CBLN2   DD_LHb   
# # 3     0       3.37 0.73  0.129         0       5 GALR1   DD_LHb  
# 
# top5_genes_habenula <- as.data.frame(top5)$gene



## ========================================================================== ##

## (1) This identifies cis-regulatory elements by linking chromatin-accessible peaks to gene expression using correlation (and optionally accounting for covariates).
## Here I am considering all the 42 WNN clusters

## ========================================================================== ##

## Find peaks that are correlated with the expression of nearby genes

## - link peaks to gene expression using correlation 
## - For each gene, LinkPeaks() computes the correlation coefficient (CC) between the gene expression and accessibility of each peak within a given distance from the gene TSS, and computes an expected CC for each peak given the GC content, accessibility, and length of the peak. The expected coefficient values for the peak are then used to compute a z-score and p-value.

message("Link Peaks to Genes (cis-regulatory analysis)")

SeuratOBJ

# Signac requires GC content or other DNA sequence information for each peak in order to compute the correlation between accessibility and expression properly. This is handled using the RegionStats() function before calling LinkPeaks().

genome <- BSgenome.Hsapiens.UCSC.hg38

SeuratOBJ <- RegionStats(
  object = SeuratOBJ,
  genome = genome,
  assay = "ATAC"  
)

SeuratOBJ[["ATAC"]]
# ChromatinAssay data with 262951 features for 55702 cells
# Variable features: 249866 
# Genome: 
# Annotation present: TRUE 
# Motifs present: FALSE 
# Fragment files: 10 

# Create pseudo-bulk replicates

granges(SeuratOBJ)
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths

# remove the features that correspond to chromosome scaffolds or other sequences instead of the (22+2) standard chromosomes
peaks.keep <- seqnames(granges(SeuratOBJ)) %in% standardChromosomes(granges(SeuratOBJ))
tryCatch(
  {
    SeuratOBJ <- SeuratOBJ[as.vector(peaks.keep), ]
  }, error = function(e) {
    message(e)
  })

# Find peaks that are correlated with the expression of nearby genes

atac <- LinkPeaks(
  object = SeuratOBJ,
  peak.assay = "ATAC",
  expression.assay = "RNA",  # Make sure RNA assay is integrated
  # genes.use = top5_genes_habenula,      # Or supply vector of gene names if you're interested in a subset
  method = "pearson",         # I am starting with default settings
  distance = 1e5             # Cis distance (e.g., 100kb window)
)

# Testing 39 genes and 262891 peaks
# Found gene coordinates for 31 genes
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=43s  

# ## Optionally find peaks by clusters and accounting for covariates
# 
# for (clus in unique(top5$cluster)) {
#   # testing: clus = 5
#   # tmp_name <- paste0(
#   #   Seurat_base_name,
#   #   "_PEAKS_hb-cluster-",
#   #   clus,
#   #   ".pdf"
#   # )
#   
#   message("Finding peaks for habenula cluster: ", clus)
#   
#   top5_genes <- top5 |>
#     filter(cluster == clus)
#   top5_genes <- as.data.frame(top5_genes)$gene
#   print(top5_genes)
#   
#   # pdf(file = here(plotDir, tmp_name))
#   # dev.off()
#   
# }
# 
# 
# # Coverage and gene link plot for a gene of interest
# CoveragePlot(
#   object = atac,
#   region = "GENE_NAME",      # Replace with e.g. "CD14"
#   features = "GENE_NAME",
#   expression.assay = "RNA",
#   extend.upstream = 10000,
#   extend.downstream = 10000
# )
                   
