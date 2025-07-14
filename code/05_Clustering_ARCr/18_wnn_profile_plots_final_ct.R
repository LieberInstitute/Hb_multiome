########################################################################
## Plot WNN with final cell types 
##
## Authors. CSC
## Date. Jun 30, 2025
##
## Recommended resources on interactive mode: srun --pty --mem=80GB --x11 bash
########################################################################

library("Seurat")
#library("SingleCellExperiment")
#library("dendextend")
#library("dynamicTreeCut")
library("ggplot2")
library("patchwork")
library("dplyr")
library("stringr")
library("here")

# directories

Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds_renamed_visium_HD.rds"
inputSeuratRDS <- here(
    "processed-data", 
    "05_Clustering_ARCr", 
    "17_wnn_clustering_final_ct", 
    Seurat_base_name)

# input_ct_summary_CSV <- here(
#     "processed-data",
#     "05_Clustering_ARCr",
#     "16_rename_idents_post_mean_ratio_HClust",
#     "WNN_full_annotation_meta_data.csv")

processedDir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "18_wnn_profile_plots_final_ct.R"
)

plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "18_wnn_profile_plots_final_ct.R"
)

## Check directories
if (!dir.exists(plotDir)) {
    dir.create(plotDir)
}
if (!dir.exists(processedDir)) {
    dir.create(processedDir)
}

#===============================================================================

message(Sys.time(), " - Load Seurat with WNN clusters")

SeuratOBJ <- readRDS(inputSeuratRDS)
SeuratOBJ
DefaultAssay(SeuratOBJ) <- "RNA"

## some verification
class(SeuratOBJ[["ATAC"]])
levels(SeuratOBJ)
colnames(SeuratOBJ@meta.data)
# check cell-types
table(SeuratOBJ$merged_cluster)
# Astrocyte       Endo Excit_Thal Inhib_Thal        LHb        MHb  Microglia 
# 2684        343       9738       6024      19673      10944        663 
# Oligo        OPC       Thal 
# 4882        638        111 

## =============================================================================


message("Processing plots on final cell-types ...")

## Some visualizations: VPlots, DimPlot, Feature, DotPlot ...

## extract suffix name to give unique name to plots
seurat_name <- str_extract(Seurat_base_name, pattern = "k[3:4]0\\_C\\.\\w*")

Reductions(SeuratOBJ)

plt1 <- DimPlot(SeuratOBJ, 
                label = TRUE, 
                reduction = "wnn.umap",
                label.size = 3) + 
    NoLegend() +
    labs(title = "WNN cell types")

ggsave(here(plotDir, "WNN_umap.pdf"), plt1, width = 7, height = 7)

plt2 <- DimPlot(SeuratOBJ, 
                label = TRUE, 
                reduction = "umap.lsi.integrated",
                label.size = 3) + 
    NoLegend() +
    labs(title = "WNN cell types in atac")

ggsave(here(plotDir, "WNN_umap_lsi_integrated.pdf"), plt2, width = 7, height = 7)

dim_plots <- (plt1 + plt2)
ggsave(here(plotDir, "Hb_DimPlot_WNN_ALL_merged_clusters_side_to_side.pdf"), dim_plots, width = 10, height = 5)


## merged clusters, I picked up Hex color codes similar to those used on human pilot
my_colors <- c(
    LHb = "#1f78b4",
    MHb = "#b74d4d",
    Oligo = "#384a08",
    Astrocyte = "#890606", 
    OPC = "#829454",
    Microglia = "#141b02",
    Endo = "#d95f02",
    Inhib_Thal = "#9a9fe7",
    Excit_Thal = "#42467b",
    Thal = "#4d55b7"
)

plt1 <- DimPlot(SeuratOBJ, 
                label = FALSE, 
                reduction = "wnn.umap",
                group.by = "merged_cluster", 
                label.size = 3,
                cols = my_colors) + 
    #NoLegend() +
    labs(title = "WNN Broad cell-types")

ggsave(here(plotDir, "WNN_merged_clusters.pdf"), plt1, width = 8, height = 7)


message("WNN UMAP done!")


# ## Compute Hierarchical Clustering on PC. ----- FASTER VERSION
# ## - means only ~10–30 dimensions, making it fast.
# ## - generate the average profiles per cluster
# ## - even randomly sample cells for a quick dendrogram
# 
# ## Pull PCA embeddings: cells x PCs
# ## Get cluster identities
# clusters <- Idents(SeuratOBJ)
# head(clusters)
# head(SeuratOBJ$cluster_ann)
# pca_mat <- Embeddings(SeuratOBJ, reduction = "pca")
# head(pca_mat)
# rm("SeuratOBJ")
# 
# ## Compute cluster centroids in PCA space
# cluster_means <- as.data.frame(pca_mat) |>
#     mutate(cluster = clusters) |>
#     group_by(cluster) |>
#     summarize(across(starts_with("PC"), mean), .groups = "drop")
# head(cluster_means)
# 
# ## Convert back to matrix (clusters x PCs)
# cluster_mat <- as.matrix(cluster_means[,-1])
# rownames(cluster_mat) <- cluster_means$cluster
# head(rownames(cluster_mat))
# 
# dend_cluster <- dist(cluster_mat) |> 
#     hclust(method = "ward.D2") |> 
#     as.dendrogram(hang = 0.2)
# 
# str(dend_cluster)
# # main branch: 'dendrogram' with 2 branches and 41 members total, at height 73.22106
# 
# ## Method: "ward.D2"
# # - improved method, mathematically consistent version of Ward’s hierarchical clustering
# # - minimizes the total within-cluster variance (the sum of squared deviations from cluster means)
# # - tells hclust to explicitly compute merges using squared Euclidean distances
# 
# ## Save data & plot
# message(Sys.time(), " - Save")
# save(dend_cluster, file = here(processedDir, "wnn_hierarchical_cluster_wnn-pca.Rdata"))
# 
# # Create categories to color branches
# cluster_means$cluster
# cluster_categories <- ifelse(grepl("LHb", cluster_means$cluster), "LHb",
#                              ifelse(grepl("MHb", cluster_means$cluster), "MHb", "No-Hb"))
# names(cluster_categories) <- cluster_means$cluster
# 
# # Map categories to colors
# category_colors <- ifelse(cluster_categories[labels(dend_cluster)] == "LHb", "tomato",
#                      ifelse(cluster_categories[labels(dend_cluster)] == "MHb", "darkblue", "black"))
# 
# message(Sys.time(), " - Plot Dendrograms - Cluster centroids in PCA")
# 
# pdf(file = here(plotDir, "dendrogram_cluster_centroid_on_pca.pdf"), width = 12, height = 8)
# 
# # Set settings 
# dend_with_heights <- dend_cluster |> 
#     set("labels_cex", 0.8) |>
#     set("labels_col", category_colors) |>
#     set("nodes_pch", 19) |>
#     set("nodes_cex", 0.7) |>
#     set("nodes_col", "blue") |>
#     set("leaves_col", "darkred") 
# 
# plot(
#     dend_cluster,
#     main = "WNN Hierarchical clustering",
#     ylab = "squared PCA distances (ward.D2)",
#     lwd = 1.5
# )
# 
# dev.off()
# 
# 
# 
# message(Sys.time(), "Dendrograms Done!")



# library("slurmjobs")
# job_single(
#     "17_wnn_hierarchical_clustering_final_ct", 
#     cores = 2, 
#     partition = "katun", 
#     memory = "80G", 
#     create_shell = TRUE
#     )

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



