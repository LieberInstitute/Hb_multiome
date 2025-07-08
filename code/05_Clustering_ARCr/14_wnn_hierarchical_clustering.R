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
head(cluster_means)

## Convert back to matrix (clusters x PCs)
cluster_mat <- as.matrix(cluster_means[,-1])
rownames(cluster_mat) <- cluster_means$cluster
head(rownames(cluster_mat))

dend_cluster <- dist(cluster_mat) |> 
    hclust(method = "ward.D2") |> 
    as.dendrogram(hang = 0.2)

str(dend_cluster)
# main branch: 'dendrogram' with 2 branches and 41 members total, at height 73.22106

## Method: "ward.D2"
# - improved method, mathematically consistent version of Ward’s hierarchical clustering
# - minimizes the total within-cluster variance (the sum of squared deviations from cluster means)
# - tells hclust to explicitly compute merges using squared Euclidean distances

## Save data & plot
message(Sys.time(), " - Save")
save(dend_cluster, file = here(processedDir, "wnn_hierarchical_cluster_wnn-pca.Rdata"))

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
    main = "WNN Hierarchical clustering",
    ylab = "squared PCA distances (ward.D2)",
    #xlab = "WNN clusters",
    cex = 0.8,
    lwd = 1.5
    )

dev.off()

#===============================================================================


## Substract pairs of leaves form all clusters to add to full summary meta-data

#str(dend_cluster)
dend_cluster %>% nnodes
# [1] 81
dend %>% nleaves
# [1] 41
dend_cluster %>% get_nodes_attr("label")
# [1] NA                  NA                  "C.21.Astrocyte"   
# [4] NA                  NA                  "C.41.Microglia"   
# [7] NA                  "C.02.Oligo"        "C.22.Oligo"   ...
dend_cluster %>% labels
# [1] "C.21.Astrocyte"    "C.41.Microglia"    "C.02.Oligo"       
# [4] "C.22.Oligo"        "C.29.Endo"         "C.27.Microglia"   
# [7] "C.20.Astrocyte"    "C.26.OPC"          "C.07.MHb.2"   

dend_cluster %>% get_nodes_attr("members")
dend_cluster %>% get_nodes_attr("members", id = c(2)) # number of members for nodes 2 and 5

# Apply to all nodes and extract all leaves at any level
all_leaves_per_node <- list()
i <- 0

extract_leaves <- function(node) {
    if (!is.leaf(node)) {
        # collect all leaf labels under the corresponding node
        leaves <- labels(node)
        if (length(leaves) > 1) {
            # store only if it's a true branch
            i <<- i + 1
            all_leaves_per_node[[i]] <<- leaves
        }
    }
    return(node)
}
dendrapply(dend_cluster, extract_leaves)

## check the character vectors with all leaves
str(all_leaves_per_node)
## only keep branches that have exactly two items (two leaves under that branch)
branches_with_two_leaves <- Filter(function(x) length(x) == 2, all_leaves_per_node)
branches_with_two_leaves

## compose a df
two_leaf_pairs_df <- lapply(
    seq_along(
        branches_with_two_leaves
        ), function(branch_id) {
    leaves <- branches_with_two_leaves[[branch_id]]
    data.frame(branch_id = branch_id,
               leaf1 = grep("^C\\.dd\\.", leaves[1]),
               leaf2 = leaves[2],
               stringsAsFactors = FALSE)
}) |> bind_rows()

two_leaf_pairs_df <- two_leaf_pairs_df |>
    mutate(
        HClust_pairs = paste(str_extract(two_leaf_pairs_df$leaf1, "\\d{2}"), "-", str_extract(two_leaf_pairs_df$leaf2, "\\d{2}"))
    )
head(two_leaf_pairs_df)
# branch_id             leaf1           leaf2 HClust_pairs
# 1          1        C.02.Oligo      C.22.Oligo      02 - 22
# 2          2    C.20.Astrocyte        C.26.OPC      20 - 26
# 3          3        C.07.MHb.2      C.36.MHb.3      07 - 36
# 4          4      C.11.MHb.1.2      C.14.MHb.1      11 - 14
# 5          5         C.13.Thal C.39.Inhib.Thal      13 - 39
# 6          6 C.06.ExcitT.LHb.4 C.04.Excit.Thal      06 - 04
# 7          7          C.24.LHb    C.30.MHb.LHb      24 - 30
# 8          8    C.18.LHb.1.3.4      C.23.LHb.1      18 - 23
# 9          9   C.28.Inhib.Thal C.38.Inhib.Thal      28 - 38
# 10        10   C.01.Inhib.Thal C.19.Inhib.Thal      01 - 19
# 11        11   C.12.Excit.Thal C.32.Excit.Thal      12 - 32
# 12        12   C.03.Excit.Thal       C.37.Thal      03 - 37
# 13        13   C.15.Excit.Thal C.35.Excit.Thal      15 - 35

## check
all_leaves_with_pairs <- sort(c(str_extract(two_leaf_pairs_df$leaf1, "\\d{2}"), str_extract(two_leaf_pairs_df$leaf2, "\\d{2}")))
all_leaves_with_pairs

## save HCLust leave-pairs
f_name <- here(processedDir, "WNN_HClust_leave_pair.csv")
write.csv(two_leaf_pairs_df, f_name, row.names = FALSE)


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



