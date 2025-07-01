########################################################################
## Compute and Plot gene expression plots for top marker genes for one cell type 
## I used Deconvobuddies::findMarkers_1vAll(), a convenient wrapped (Scran/Deconvobuddies) to compute test.type="binom" (1vsALL)
## - Used: Default direction = "up".  Impacts p-values: if "up" genes with logFC < 0 will have p.value = 1
##
## Authors. CSC
## Date. Jun 30, 2025
##
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
## DeconvoBuddies 1.1+ is not available for Bioconductor 3.18 (conda_R/4.3.x), so I am using conda_R/4.4.x instead
########################################################################

library("SingleCellExperiment")
library("dplyr")
library("ggplot2")
library("here")


## directories
inputSCE_Dir <- here(
    "processed-data", 
    "08_spatial_registration_vs_multiome_snRNA-seq"
)
processedDir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "14_wnn_hierarchical_clustering"
)
plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "14_wnn_hierarchical_clustering"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(processedDir)) {
    dir.create(processedDir)
}

#===============================================================================

## SingleCellExperiment object derived from Seurat rna modality
sce_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v4.rds"
sce_name <- here(inputSCE_Dir, sce_base_name)
sce_name
title_name <- str_extract(sce_base_name, regex("C\\.\\w*\\_r2"))
title_name

## Load SCE derived from Seurat WNN
sce <- readRDS(sce_name)

message("====== Verifications ======")
sce
# class: SingleCellExperiment 
# dim: 29690 55702 
# ...

## Verification
head(rownames(sce))
colnames(colData(sce))

## check cell types and counts
cluster_counts_df <- as.data.frame(table(sce$cluster_ann))
head(cluster_counts_df)

# check rowData ensembl names and gene id(s)
head(rowData(sce)$gene_symbol)
head(rowData(sce)$gene_id)

celltype_counts <- table(colData(sce)$cluster_ann)
as.data.frame(celltype_counts)


#===============================================================================

library("dendextend")
library("dynamicTreeCut")


## Perform hierarchical clustering based on relative expression patterns across genes for exploratory single-cell analysis

# Compuse dist.mtx and scale each gene to make clustering based on patterns of up/down regulation, not on raw magnitude [range 0-1]
dend <- dist(t(logcounts(sce))) |> 
    scale() |> 
    hclust(mat_scaled, "ward.D2") |>
    as.dendrogram(hang = 0.2) 

dend |> unclass |> str

dend |> head

dend |> plot


#===============================================================================
## Compute Hierarchical Clustering for both logcounts (for magnitude differences) and scale data (highlights relative patterns)

message(Sys.time(), " - Cluster Dendrogram from logcounts (no scaling) ")

mat_logcounts <- t(logcounts(sce)) # cells as rows, genes as columns
dist_logcounts <- dist(mat_logcounts)
hc_logcounts <- hclust(dist_logcounts, method = "ward.D2")
dend_logcounts <- as.dendrogram(hc_logcounts, hang = 0.2)
dend_logcounts <- set(dend_logcounts, "branches_col", "blue") # optional color for clarity

## Save data
message(Sys.time(), " - Save")
save(dend_logcounts, file = here(processedDir, "wnn_hierarchical_cluster_logcounts.Rdata"))



message(Sys.time(), " - Cluster Dendrogram from scaled data ")

mat_scaled <- scale(mat_logcounts) # scale genes to mean=0, sd=1
dist_scaled <- dist(mat_scaled)
hc_scaled <- hclust(dist_scaled, method = "ward.D2")
dend_scaled <- as.dendrogram(hc_scaled, hang = 0.2)
dend_scaled <- set(dend_scaled, "branches_col", "red")

## Save data
message(Sys.time(), " - Save")
save(dend_scaled, file = here(processedDir, "wnn_hierarchical_cluster_scale_data.Rdata"))

#===============================================================================

message(Sys.time(), " - Plot Dendrograms Side by side comparison")

pdf(file = here("plots", "dendrogram_comparison_logcounts_vs_scaled.pdf"), width = 12, height = 8)

tanglegram(dend_logcounts, dend_scaled,
           main_left = "Logcounts (unscaled)",
           main_right = "Scaled data (gene z-scores)",
           common_subtrees_color_lines = TRUE,
           highlight_distinct_edges = TRUE,
           columns_width = c(5,5),
           lab.cex = 0.5)

dev.off()


# library("slurmjobs")
# job_single("14_wnn_hierarchical_clustering", cores = 2, partition = "katun", create_shell = TRUE)





