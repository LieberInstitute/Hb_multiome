#Running scdblfinder to check where called doublets fall on umap



library(SingleCellExperiment)
library(scDblFinder)
library(Seurat)
library(MetaMarkers)
library(dplyr)
library(ggplot2)
library(scater)
library(scran)
library(here)

here::here()


multiome_path = here('processed-data', '08_spatial_registration_vs_multiome_snRNA-seq','mid')

#Multiome human data
multiome_sce = readRDS(paste0(multiome_path, '/seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v5.rds'))
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))

colnames(colData(multiome_sce))

table(multiome_sce$orig.ident)
#Run scDblFinder defaults

multiome_sce <- scDblFinder(multiome_sce, samples = multiome_sce$orig.ident)

table(multiome_sce$scDblFinder.class)
#singlet doublet 
#  50208    5308 


multiome_sce_meta = colData(multiome_sce)

#This is the donor integrated LHb4 and 7, has the putative inhibitory cluster annotations
script_02_data_path = here('processed-data', '99_donor_cluster_replicability','02_LHb4_investigation')
multiome_seurat_integrated = readRDS(paste0(script_02_data_path, '/multiome_LHb4_LHb7_integrated_seurat.rds'))

inhib_meta = multiome_seurat_integrated[[]]

# Add the putative Inhibitory neuron annotations to the midSeurat object
inhib_barcodes_1 = rownames(inhib_meta)[inhib_meta$refined_mid_cluster == 'Putative_Inhib_LHb_1']
inhib_barcodes_2 = rownames(inhib_meta)[inhib_meta$refined_mid_cluster == 'Putative_Inhib_LHb_2']

# Add labels
multiome_sce$refined_mid_cluster = multiome_sce$mid_cluster
multiome_sce$refined_mid_cluster[rownames(multiome_sce_meta) %in% inhib_barcodes_1] = 'Putative_Inhib_LHb_4.1'
multiome_sce$refined_mid_cluster[rownames(multiome_sce_meta) %in% inhib_barcodes_2] = 'Putative_Inhib_LHb_4.2'

table(multiome_sce$refined_mid_cluster, multiome_sce$scDblFinder.class)
  #                      singlet doublet
  # Astrocyte               2305     379
  # Endo                     306      37
  # Excit.Thal              8706    1032
  # Inhib.Thal              5507     517
  # LHb.1                   1160     109
  # LHb.1.3                  179       7
  # LHb.1.3.4               1425     110
  # LHb.2.7                 2570     201
  # LHb.4                  10111    1593
  # LHb.7                     81     130
  # MHb.1                   4160     210
  # MHb.1.2                 3467     352
  # MHb.2                   2471     139
  # MHb.3                    123      22
  # Microglia                646      17
  # Oligo                   4553     145
  # OPC                      616      22
  # Putative_Inhib_LHb_1    1011     115
  # Putative_Inhib_LHb_2     807      64
  # Thal                       4     107



#Visualize doublets on the original umap

## Seurat object with the mid-level cluster annots and dim reductions saved
midSeurat_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "22_add_mid_level_clustering"
)
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"

mid_file_name <- here(midSeurat_Dir, Seurat_base_name)

midSeurat = readRDS(mid_file_name)

#Cell barcodes match
table(colnames(midSeurat) == rownames(colData(multiome_sce)))

#Add the sclDoubFinder class to the midSeurat object
midSeurat$scDblFinder.class = multiome_sce$scDblFinder.class
midSeurat$scDblFinder.score = multiome_sce$scDblFinder.score
midSeurat$refined_mid_cluster = multiome_sce$refined_mid_cluster

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
cluster_levels
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

message("Processing UMAP ...")

## extract suffix name to give unique name to plots
seurat_name <- stringr::str_extract(Seurat_base_name, pattern = "k[3:4]0\\_C\\.\\w*")


#Verify original UMAP
plt1 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "mid_cluster",
                label.size = 3,
                cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "WNN cell types (Mid-resolution)")

#Looks good
plt1

my_colors_mid
my_colors_mid["Putative_Inhib_LHb_4.1"] <- "#8B0000"  # Dark red
my_colors_mid["Putative_Inhib_LHb_4.2"] <- "#DC143C"  # Crimson red


#Verify original UMAP
plt2 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "refined_mid_cluster",
                label.size = 3,
                cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "WNN cell types (Mid-resolution)")

plt2

#Verify original UMAP
plt3 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "scDblFinder.class",
                label.size = 3) + 
    labs(title = "WNN cell types (Mid-resolution)")

plt3

FeaturePlot(midSeurat, features = "scDblFinder.score", reduction = "wnn.umap")

cluster_order <- midSeurat@meta.data |>
  dplyr::count(refined_mid_cluster, scDblFinder.class) |>
  dplyr::group_by(refined_mid_cluster) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::filter(scDblFinder.class == "doublet") |>
  dplyr::arrange(prop) |>
  dplyr::pull(refined_mid_cluster)

midSeurat@meta.data |>
  dplyr::mutate(refined_mid_cluster = factor(refined_mid_cluster, levels = cluster_order)) |>
  ggplot(aes(x = refined_mid_cluster, fill = scDblFinder.class)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Refined mid cluster", y = "Proportion", fill = "scDblFinder") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



#Andquick check in the integrated LHb4 and Lhb7 data

multiome_seurat_integrated

index = match(colnames(multiome_seurat_integrated), colnames(multiome_sce))

multiome_seurat_integrated$sclDblFinder.class = multiome_sce$scDblFinder.class[index]
multiome_seurat_integrated$sclDblFinder.score = multiome_sce$scDblFinder.score[index]

DimPlot(multiome_seurat_integrated, group.by = 'sclDblFinder.class', reduction = 'umap', pt.size = 1)
DimPlot(multiome_seurat_integrated, group.by = 'seurat_clusters', reduction = 'umap', pt.size = 1)
DimPlot(multiome_seurat_integrated, group.by = 'refined_mid_cluster', reduction = 'umap', pt.size = 1)

FeaturePlot(multiome_seurat_integrated, features = "GPR151", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'CPM')
FeaturePlot(multiome_seurat_integrated, features = "POU4F1", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'CPM')



