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

library("Seurat")
library("SingleCellExperiment")
library("dendextend")
library("dynamicTreeCut")
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
# title_name <- str_extract(sce_base_name, regex("C\\.\\w*\\_r2"))
# title_name

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

## Compute Hierarchical Clustering on PC. ----- FASTER VERSION
## - means only ~10–30 dimensions, making it fast.
## - generate the average profiles per cluster
## - even randomly sample cells for a quick dendrogram

message(Sys.time(), " - Cluster Dendrogram on WNN clusters on PCA space")

outputRDS_Dir <- here("processed-data", "05_Clustering_ARCr", "05_rename_idents")
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2"
Seurat_base_name <- paste0(Seurat_base_name, "_renamed_visium.rds")
# Load Seurat with WNN idents given by default 
seurat_RDSname <- here(outputRDS_Dir, Seurat_base_name)
SeuratOBJ <- readRDS(seurat_RDSname)
SeuratOBJ
DefaultAssay(SeuratOBJ) <- "RNA"

## Pull PCA embeddings: cells x PCs
## Get cluster identities
clusters <- Idents(SeuratOBJ)
head(clusters)
head(SeuratOBJ$cluster_ann)
pca_mat <- Embeddings(SeuratOBJ, reduction = "pca")
head(pca_mat)
rm("SeuratOBJ")

## Compute cluster centroids in PCA space
cluster_means <- as.data.frame(pca_mat) |>
    mutate(cluster = clusters) |>
    group_by(cluster) |>
    summarize(across(starts_with("PC"), mean), .groups = "drop")

# Convert back to matrix (clusters x PCs)
cluster_mat <- as.matrix(cluster_means[,-1])
rownames(cluster_mat) <- cluster_means$cluster
head(rownames(cluster_mat))

dist_cluster <- dist(cluster_mat)
hc_cluster <- hclust(dist_cluster, method = "ward.D2")
dend_cluster <- as.dendrogram(hc_cluster, hang = 0.2)

## Save data & plot
message(Sys.time(), " - Save")
#save(dend_pca, file = here(processedDir, "wnn_hierarchical_cluster_rna-pca.Rdata"))

# Create categories to color branches
cluster_means$cluster
cluster_categories <- ifelse(grepl("LHb", cluster_means$cluster), "LHb",
                             ifelse(grepl("MHb", cluster_means$cluster), "MHb", "No-Hb"))
names(cluster_categories) <- cluster_means$cluster

# Map categories to colors
category_colors <- ifelse(cluster_categories[labels(dend_cluster)] == "LHb", "tomato",
                     ifelse(cluster_categories[labels(dend_cluster)] == "MHb", "darkblue", "black"))

# Set colors
dend_cluster <- dend_cluster |> set("labels_col", category_colors)


message(Sys.time(), " - Plot Dendrograms - Cluster centroids in PCA")

pdf(file = here(plotDir, "dendrogram_cluster_centroid_on_pca.pdf"), width = 12, height = 8)

plot(
    dend_cluster, 
    main = "Hierarchical clustering of WNN cluster centroids (PCA space)",
    ylab = "Height",
    cex = 0.8,
    lwd = 1.5
    )

dev.off()


#===============================================================================

## Compute Hierarchical Clustering for both logcounts (for magnitude differences) and scale data (highlights relative patterns)

message(Sys.time(), " - Cluster Dendrogram from logcounts (no scaling) ")

mat_logcounts <- t(logcounts(sce)) # cells as rows, genes as columns
dist_logcounts <- dist(mat_logcounts)
hc_logcounts <- hclust(dist_logcounts, method = "ward.D2")
dend_logcounts <- as.dendrogram(hc_logcounts, hang = 0.2)
dend_logcounts <- set(dend_logcounts, "branches_col", "blue") # optional color for clarity

## Save data
#message(Sys.time(), " - Save")
#save(dend_logcounts, file = here(processedDir, "wnn_hierarchical_cluster_logcounts.Rdata"))

pdf(file = here(plotDir, "dendrogram_wnn_logcount.pdf"), width = 12, height = 8)

plot(
    dend_logcounts, 
    main = "Hierarchical clustering of logcounts",
    ylab = "Height",
    cex = 0.8,
    lwd = 1.5
)

dev.off()


message(Sys.time(), " - Cluster Dendrogram from scaled data ")

mat_scaled <- scale(mat_logcounts) # scale genes to mean=0, sd=1
dist_scaled <- dist(mat_scaled)
hc_scaled <- hclust(dist_scaled, method = "ward.D2")
dend_scaled <- as.dendrogram(hc_scaled, hang = 0.2)
dend_scaled <- set(dend_scaled, "branches_col", "red")

## Save data
#message(Sys.time(), " - Save")
#save(dend_scaled, file = here(processedDir, "wnn_hierarchical_cluster_scale_data.Rdata"))

pdf(file = here(plotDir, "dendrogram_wnn_scaledata.pdf"), width = 12, height = 8)

plot(
    dend_scaled, 
    main = "Hierarchical clustering of scale-data",
    ylab = "Height",
    cex = 0.8,
    lwd = 1.5
)

dev.off()

 
# #===============================================================================
 
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


message(Sys.time(), "Dendrograms Done!")



# library("slurmjobs")
# job_single("14_wnn_hierarchical_clustering", cores = 2, partition = "katun", create_shell = TRUE)

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



