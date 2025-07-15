########################################################################
## Plot WNN with final cell types 
##
## Authors. CSC
## Date. Jun 30, 2025
##
## Recommended resources on interactive mode: srun --pty --mem=80GB --x11 bash
########################################################################

library("Seurat")
library("dendextend")
library("dynamicTreeCut")
library("purrr")
library("ggplot2")
library("ggtext") # Build names with HTML color tags / DotPlot
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
table(SeuratOBJ$seurat_clusters)
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 3906 3371 3271 2911 2771 2667 2610 2537 2430 2344 2212 2187 2036 2026 1625 1607 
# 17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32 
# 1587 1535 1380 1343 1341 1327 1269  825  707  638  587  543  343  213  209  196 
# 33   34   35   36   37   38   39   40   41 
# 186  184  165  145  111  105   90   84   76 

table(sort(SeuratOBJ$cluster_ann))
# C.01.Inhib.Thal        C.02.Oligo   C.03.Excit.Thal   C.04.Excit.Thal 
# 3906              3371              3271              2911 
# C.05.LHb.2.7 C.06.ExcitT.LHb.4        C.07.MHb.2        C.08.LHb.4 
# 2771              2667              2610              2537 

table(SeuratOBJ$merged_cluster)
# Astrocyte       Endo Excit_Thal Inhib_Thal        LHb        MHb  Microglia 
# 2684        343       9738       6024      19673      10944        663 
# Oligo        OPC       Thal 
# 4882        638        111 

levels(SeuratOBJ)
# [1] "C.04.LHb.4"      "C.05.LHb.2.7"    "C.06.LHb.4"      "C.07.MHb.2"     
# [5] "C.08.LHb.4"      "C.09.LHb.4"      "C.10.MHb.1"      "C.11.MHb.1.2"   
# [9] "C.13.LHb.4"      "C.14.MHb.1"      "C.16.MHb.1.2"    "C.18.LHb.1.3.4" 
# [13] "C.23.LHb.1"      "C.24.LHb.4"      "C.30.LHb.7"      "C.31.LHb.4"     
# [17] "C.33.LHb.1.3"    "C.36.MHb.3"      "C.40.LHb.4"      "C.01.Inhib.Thal" ....

## =============================================================================


message("Processing UMAP plots ...")

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


plt1 <- DimPlot(SeuratOBJ, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "seurat_clusters",
                label.size = 3) + 
    NoLegend() +
    labs(title = "WNN cell types")

ggsave(here(plotDir, "WNN_umap_clusterID.pdf"), plt1, width = 7, height = 7)


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
    MHb = "#ad1d8c",
    Oligo = "#384a08",
    Astrocyte = "#532222", 
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

## =============================================================================


message("Processing Violin plots ...")

## Plot Hb canonical genes for merged_clusters

plot_violin_merged_clusters <- function(seurat_obj, genes, group_col = "merged_cluster", colors = NULL) {
    
    plots <- purrr::map(genes, ~ {
        VlnPlot(
            object = SeuratOBJ,
            slot = "data",
            group.by = group_col, 
            features = .x,
            pt.size = 0.2,
            alpha = 0.1,
            cols = my_colors
        ) +
            labs(title = .x) +
            theme(
                text = element_text(size = 10),
                axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1),
                axis.text.y = element_text(size = 10),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 12)
            ) +
            coord_flip() +
            NoLegend()
    }) |> purrr::set_names(genes)
    
    return(plots)
    
}

genes_to_plot <- c("GPR151", "POU4F1", "TAC3")
my_plots <- plot_violin_merged_clusters(SeuratOBJ, genes_to_plot, colors = my_colors)

plt1 <- my_plots[["GPR151"]] + my_plots[["POU4F1"]] + my_plots[["TAC3"]] 

ggsave(here(plotDir, "WNN_Vplots_Hb_canonical_merged_clusters.pdf"), plt1, width = 6, height = 7)

message("WNN UMAP done!")

## =============================================================================


message("Processing GeneExpression Dot plots ...")

# Build color mapping: LHb and MHb get colors, others default to black
label_colors <- ifelse(grepl("LHb", clusters), "#1f78b4", 
                       ifelse(grepl("MHb", clusters), "#ad1d8c",
                              "black"))
# Build names with HTML color tags
clusters_colored <- paste0("<span style='color:", label_colors, "'>", clusters, "</span>")
names(clusters_colored) <- clusters  # keep mapping

plt1 <- DotPlot(SeuratOBJ, 
        features = genes_to_plot) +
        #group.by = "merged_cluster") +
    theme(
        text = element_text(size = 12),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_markdown(size = 9),  # ggtext to parse html
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
    )  +
    scale_y_discrete(labels = clusters_colored)  # apply colored labels

ggsave(here(plotDir, "WNN_DotPlot_Hb_canonical_all_clusters.pdf"), plt1, width = 5, height = 7)




## =============================================================================

message("Processing HClust plots ...")

## Compute Hierarchical Clustering on PC. ----- FASTER VERSION
## - means only ~10–30 dimensions, making it fast.
## - generate the average profiles per cluster
## - even randomly sample cells for a quick dendrogram

## Pull PCA embeddings: cells x PCs
## Get cluster identities
clusters <- Idents(SeuratOBJ)
table(clusters)
head(clusters)
head(SeuratOBJ$cluster_ann)
pca_mat <- Embeddings(SeuratOBJ, reduction = "pca")
head(pca_mat)
# rm("SeuratOBJ")

## Compute cluster centroids in PCA space
cluster_means <- as.data.frame(pca_mat) |>
    mutate(cluster = clusters) |>
    group_by(cluster) |>
    summarize(across(starts_with("PC"), mean), .groups = "drop")
head(cluster_means)
 
## Convert back to matrix (clusters x PCs)
cluster_mat <- as.matrix(cluster_means[,-1])
rownames(cluster_mat) <- cluster_means$cluster
head(rownames(cluster_mat))

dend_cluster <- dist(cluster_mat) |>
    hclust(method = "ward.D2") |>
    as.dendrogram(hang = 0.2)

#str(dend_cluster)
# main branch: 'dendrogram' with 2 branches and 41 members total, at height 73.22106
 
## Method: "ward.D2"
# - improved method, mathematically consistent version of Ward’s hierarchical clustering
# - minimizes the total within-cluster variance (the sum of squared deviations from cluster means)
# - tells hclust to explicitly compute merges using squared Euclidean distances
 
## Save data & plot
message(Sys.time(), " - Save")
#save(dend_cluster, file = here(processedDir, "wnn_hierarchical_cluster_wnn-pca.Rdata"))

# Create categories to color branches based on 
head(cluster_means$cluster)

## merged clusters, I picked up Hex color codes similar to those used on human pilot
cluster_means <- cluster_means |>
    mutate(category = case_when(
        grepl("LHb", cluster) ~ "LHb",
        grepl("MHb", cluster) ~ "MHb",
        grepl("Oligo", cluster) ~ "Oligo",
        grepl("Astrocyte", cluster) ~ "Astrocyte",
        grepl("OPC", cluster) ~ "OPC",
        grepl("Microglia", cluster) ~ "Microglia",
        grepl("Endo", cluster) ~ "Endo",
        grepl("Inhib", cluster) ~ "Inhib_Thal",
        grepl("Excit", cluster) ~ "Excit_Thal",
        grepl("Thal", cluster) ~ "Thal",
        TRUE ~ "Other"
    ))


# Named vector mapping each cluster to its category
cluster_categories <- cluster_means$category
names(cluster_categories) <- cluster_means$cluster

# Map labels on the dendrogram to categories, then to hex colors
label_categories <- cluster_categories[labels(dend_cluster)]
label_colors <- my_colors[label_categories]

message(Sys.time(), " - Plot Dendrogram - Cluster centroids in PCA")

pdf(file = here(plotDir, "dendrogram_cluster_centroid_on_pca.pdf"), width = 10, height = 4)
 
# Set settings and color vector
dend_with_heights <- dend_cluster |>
    set("labels_cex", 0.8) |>
    set("labels_col", label_colors) |>
    set("nodes_pch", 19) |>
    set("nodes_cex", 0.7) # |>
    # set("nodes_col", "blue") |>
    # set("leaves_col", "darkred")

# Rotate the tree to change orientation
dend_flipped <- rotate(dend_with_heights, order = rev(labels(dend_with_heights)))

plot(
    dend_with_heights,
    #dend_flipped,
    #horiz = TRUE,
    #main = "WNN Hierarchical clustering",
    ylab = "squared PCA distances (ward.D2)",
    lwd = 1.5
) +
    theme(plot.margin = margin(10, 10, 30, 10))  # t, r, b, l


dev.off()

message(Sys.time(), "Dendrograms Done!")



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



