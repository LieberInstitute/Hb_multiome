
#Computing doublet scores and checking the doublet proportions of each fine resolution cluster
#Drop any doublet enriched clusters
#Redo the dimensionality reduction
#Save for downstream analyses

library(SingleCellExperiment)
library(scDblFinder)
library(Seurat)
library(dplyr)
library(ggplot2)
library(here)

here::here()

#Path for new data generated
new_data_path = here('processed-data', '05_5_drop_doublets','01_drop_doublets_and_reDimReduce')
#Path to plot directory
plot_path = here('plots', '05_5_drop_doublets', '01_drop_doublets_and_reDimReduce')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


## Seurat object with the mid-level cluster annots and dim reductions saved
midSeurat_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "22_add_mid_level_clustering"
)
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"

mid_file_name <- here(midSeurat_Dir, Seurat_base_name)

midSeurat = readRDS(mid_file_name)


#Run scDblFinder to compute doublet scores
sce <- scDblFinder(GetAssayData(midSeurat , slot="counts"), samples = midSeurat$orig.ident)
# port the resulting scores back to the Seurat object:
midSeurat$scDblFinder.score <- sce$scDblFinder.score
midSeurat$scDblFinder.class <- sce$scDblFinder.class

table(midSeurat$scDblFinder.class)
#singlet doublet 
#  50184    5332 

#Fine resolution clusters
cluster_order <- midSeurat@meta.data |>
  dplyr::count(cluster_ann, scDblFinder.class) |>
  dplyr::group_by(cluster_ann) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::filter(scDblFinder.class == "doublet") |>
  dplyr::arrange(prop) |>
  dplyr::pull(cluster_ann)

midSeurat@meta.data |>
  dplyr::mutate(cluster_ann = factor(cluster_ann, levels = cluster_order)) |>
  ggplot(aes(x = cluster_ann, fill = scDblFinder.class)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Fine clusters", y = "Proportion", fill = "scDblFinder") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

midSeurat@meta.data |>
  dplyr::mutate(cluster_ann = factor(cluster_ann, levels = cluster_order)) |>
  ggplot(aes(x = cluster_ann, y = scDblFinder.score)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Fine clusters") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


#Check out current UMAP and the clusters we will drop
## =============================================================================
## Picked up Hex-color codes similar across cell-type

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

## assign color gradients to mid resolution clusters based on Broad cell-types

# extract LHb and MHb clusters
cluster_levels <- levels(midSeurat)
LHb_clusters <- grep("LHb", cluster_levels, value = TRUE)
MHb_clusters <- grep("MHb", cluster_levels, value = TRUE)

# Create tonal gradients for LHb and MHb
LHb_colors <- colorspace::sequential_hcl(length(LHb_clusters), h = 210, c = 80, l = c(30, 80))
MHb_colors <- colorspace::sequential_hcl(length(MHb_clusters), h = 320, c = 80, l = c(30, 80))

# Build full cluster color map
my_colors_mid <- setNames(rep("#bdbdbd", length(cluster_levels)), cluster_levels)
my_colors_mid[LHb_clusters] <- LHb_colors
my_colors_mid[MHb_clusters] <- MHb_colors

# assign base color for other types from your existing palette
for (category in c("Oligo", "Astrocyte", "OPC", "Microglia", "Endo", "Inhib.Thal", "Excit.Thal", "Thal")) {
    matched <- grep(category, cluster_levels, value = TRUE)
    my_colors_mid[matched] <- my_colors[[gsub("\\.", "_", category)]]
}

## =============================================================================

plt1 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "mid_cluster",
                label.size = 3,
                cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "WNN cell types (Mid-resolution)")

plt1

#Based on cluster proportions of doublets, these are the ones we would exclude
#c('C.37.Thal','C.39.Inhib.Thal','C.24.LHb.4','C.30.LHb.7','C.40.LHb.4')
midSeurat$doublet_exclude = midSeurat$mid_cluster
midSeurat$doublet_exclude[midSeurat$cluster_ann %in% c('C.37.Thal','C.39.Inhib.Thal','C.24.LHb.4','C.30.LHb.7','C.40.LHb.4')] = 'Exclude'

my_colors_mid["Exclude"] <- "#ff3b3bff"

#Verify original UMAP
plt2 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "doublet_exclude",
                label.size = 3,
                cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "WNN cell types (Mid-resolution)")


plt2


ggsave(plt1, filename = 'original_mid_res_WNN_umap.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)
ggsave(plt2, filename = 'doubletRed_mid_res_WNN_umap.pdf', path = plot_path, device = 'pdf',
width = 6, height = 5, useDingbats = FALSE)


#Drop the doublet enriched clusters
midSeurat = subset(midSeurat, doublet_exclude != 'Exclude')

table(midSeurat$cluster_ann)


#Redo essential steps, normalization, variable feature selection, scaling, integrating, and dimensionality reductions
#Just start from counts again to simplify, drop the existing dimensionality reductions.

drop_cols <- c("ATAC.weight", "RNA.weight", "doublet_exclude", "RNA_snn_res.1", 
"C.leiden", "C.leiden_atac", "C.leiden_wnn", "unintegrated_clusters" ,  "seurat_clusters" )

# Keep everything else
midSeurat@meta.data <- midSeurat@meta.data |>
  dplyr::select(-dplyr::any_of(drop_cols))


colnames(midSeurat[[]])

#Drop all the current reductions
# See current reductions
Reductions(midSeurat)

# Drop all dimensional reductions
midSeurat@reductions <- list()

# Confirm
Reductions(midSeurat)


#Rescale the data 






