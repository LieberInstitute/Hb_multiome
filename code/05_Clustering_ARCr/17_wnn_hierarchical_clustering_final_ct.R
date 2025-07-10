########################################################################
## Compute and Plot gene expression plots for top marker genes for one cell type 
##
## Authors. CSC
## Date. Jun 30, 2025
##
## Recommended resources on interactive mode: srun --pty --mem=60GB --x11 bash
########################################################################

library("Seurat")
library("SingleCellExperiment")
library("dendextend")
library("dynamicTreeCut")
library("dplyr")
library("ggplot2")
library("stringr")
library("here")

# directories
# inputSCE_Dir <- here(
#     "processed-data",
#     "08_spatial_registration_vs_multiome_snRNA-seq"
# )
# processedDir <- here(
#     "processed-data",
#     "05_Clustering_ARCr",
#     "14_wnn_hierarchical_clustering"
# )
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium.rds"
inputSeuratRDS <- here(
    "processed-data", 
    "05_Clustering_ARCr", 
    "05_rename_idents", 
    Seurat_base_name)

input_ct_summary_CSV <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "16_rename_idents_post_mean_ratio_HClust",
    "WNN_full_annotation_meta_data.csv")

plotDir <- here(
    "plots",
    "05_Clustering_ARCr",
    "17_wnn_hierarchical_clustering_final_ct"
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
current_ct <- levels(SeuratOBJ)
current_ct <- sort(current_ct)

#===============================================================================


message(Sys.time(), " - Load full annotation to renmmae Seurat idents with final cell-types")

summary_ct_df <- read.csv(input_ct_summary_CSV)
colnames(summary_ct_df)

## verify clustering order and ct assignation
summary_ct_df_sorted <- summary_ct_df |>
    arrange(cluster) |>
    select(cluster, ct_final)
summary_ct_df_sorted
## remove (*) and compose and update ct    
new_ct <- paste0(summary_ct_df$cluster, ".", summary_ct_df$ct_final)
new_ct <- gsub("\\*", "", new_ct)

## verify
cell_types_df <- data.frame(
    old_ct = current_ct,
    new_ct = new_ct
)
cell_types_df
# > cell_types_df
#           old_ct          new_ct
# 1    C.01.Inhib.Thal C.01.Inhib.Thal
# 2         C.02.Oligo      C.02.Oligo
# 3    C.03.Excit.Thal C.03.Excit.Thal
# 4    C.04.Excit.Thal      C.04.LHb.4
# 5       C.05.LHb.2.7    C.05.LHb.2.7
# 6  C.06.ExcitT.LHb.4      C.06.LHb.4
# 7         C.07.MHb.2      C.07.MHb.2
# 8         C.08.LHb.4      C.08.LHb.4
# 9  C.09.ExcitT.LHb.4      C.09.LHb.4


names(new_ct) <- current_ct
## rename idents with new cell-types
SeuratOBJ <- RenameIdents(SeuratOBJ, new_ct)
levels(SeuratOBJ)
#head(Idents(SeuratOBJ))


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
# job_single("14_wnn_hierarchical_clustering", cores = 2, partition = "katun", create_shell = TRUE)

## Reproducibility information
library("sessioninfo")
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()



