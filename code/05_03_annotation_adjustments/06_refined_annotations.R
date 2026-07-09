#Consolidating the annotation changes and saving, this will be used for all downstream analyses

#Annotation changes
#C.11.MHb.1.2 changes from MHb.1.2 to MHb.2
#LHb.4 gets split into LHb.4 and inhib.LHb.4.1 and inhib.LHb.4.2
#LHb.1, LHb.1.3, and LHb.1.3.4 all get merged into LHb.1.3.4
#C.21.Astrocyte changes to Ependymal

#We're also moving from number labels to letters
#MHb.1 - MHb_A
#MHb.2 - MHb_B
#MHb1.2 - MHb_C
#MHb.3 - MHb_D
#LHb.2.7 - LHb_A
#LHb.1.3.4 - LHb_B
#LHb.4 - LHb_C

#Make sure the mid resolution annotation changes are reflected in the cluster_ann (fine resolution) annotations


library(SingleCellExperiment)
library(Seurat)
library(qs2)
library(dplyr)
library(ggplot2)
library(here)


here::here()

#Path for new data generated
new_data_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
#Path to plot directory
plot_path = here('plots', '05_03_annotation_adjustments','06_refined_annotations')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '05_01_drop_doublets','01_drop_doublets_and_reDimReduce')

#Multiome human data, singlecellexperiment
multiome_sce = qs_read(paste0(multiome_path, '/reprocessed_doubletRemoved_multiomeHab_SCE.qs2'))

## Seurat object 
midSeurat = qs_read(paste0(multiome_path, '/reprocessed_doubletRemoved_multiomeHab_seurat.qs2'))

#This has the inhibitory LHb annotations
script_02_data_path = here('processed-data', '05_03_annotation_adjustments','02_LHb4_investigation')
multiome_seurat_integrated = readRDS(paste0(script_02_data_path, '/multiome_LHb4_LHb7_integrated_seurat.rds'))
colnames(multiome_seurat_integrated[[]])

inhib_meta = multiome_seurat_integrated[[]]

#Starting point
midSeurat$refined_mid_cluster = midSeurat$mid_cluster
midSeurat$refined_cluster_ann = as.character(midSeurat$cluster_ann)

multiome_sce$refined_mid_cluster = multiome_sce$mid_cluster
multiome_sce$refined_cluster_ann = as.character(multiome_sce$cluster_ann)

#First merge the LHb.1, LHb.1.3, and LHb.1.3.4 clusters into one cluster called LHb.1.3.4
midSeurat$refined_mid_cluster[midSeurat$refined_mid_cluster %in% c('LHb.1', 'LHb.1.3')] = 'LHb.1.3.4'
#table(midSeurat$refined_mid_cluster, midSeurat$mid_cluster)

multiome_sce$refined_mid_cluster[multiome_sce$refined_mid_cluster %in% c('LHb.1', 'LHb.1.3')] = 'LHb.1.3.4'
#table(multiome_sce$refined_mid_cluster,multiome_sce$mid_cluster)

#Now change the MHb.1.2 cluster to MHb.2
midSeurat$refined_mid_cluster[midSeurat$cluster_ann == 'C.11.MHb.1.2'] = 'MHb.2'
#table(midSeurat$refined_mid_cluster, midSeurat$mid_cluster)

multiome_sce$refined_mid_cluster[multiome_sce$cluster_ann == 'C.11.MHb.1.2'] = 'MHb.2'
#table(multiome_sce$refined_mid_cluster, multiome_sce$mid_cluster)


#And update the C.21.Astrocyte to Ependymal
#Fine cluster first
midSeurat$refined_cluster_ann[midSeurat$cluster_ann == 'C.21.Astrocyte'] = 'C.21.Ependymal'
multiome_sce$refined_cluster_ann[multiome_sce$cluster_ann == 'C.21.Astrocyte'] = 'C.21.Ependymal'

#And match at the mid resolution
midSeurat$refined_mid_cluster[midSeurat$cluster_ann == 'C.21.Astrocyte'] = 'Ependymal'
multiome_sce$refined_mid_cluster[multiome_sce$cluster_ann == 'C.21.Astrocyte'] = 'Ependymal'


###############
#Now switching to letter labels for the mid resolution clusters
#
###############
midSeurat$refined_mid_cluster[midSeurat$refined_mid_cluster == 'MHb.1'] = 'MHb_A'
midSeurat$refined_mid_cluster[midSeurat$refined_mid_cluster == 'MHb.2'] = 'MHb_B'
midSeurat$refined_mid_cluster[midSeurat$refined_mid_cluster == 'MHb.1.2'] = 'MHb_C'
midSeurat$refined_mid_cluster[midSeurat$refined_mid_cluster == 'MHb.3'] = 'MHb_D'
midSeurat$refined_mid_cluster[midSeurat$refined_mid_cluster == 'LHb.2.7'] = 'LHb_A'
midSeurat$refined_mid_cluster[midSeurat$refined_mid_cluster == 'LHb.1.3.4'] = 'LHb_B'
midSeurat$refined_mid_cluster[midSeurat$refined_mid_cluster == 'LHb.4'] = 'LHb_C'

multiome_sce$refined_mid_cluster[multiome_sce$refined_mid_cluster == 'MHb.1'] = 'MHb_A'
multiome_sce$refined_mid_cluster[multiome_sce$refined_mid_cluster == 'MHb.2'] = 'MHb_B'
multiome_sce$refined_mid_cluster[multiome_sce$refined_mid_cluster == 'MHb.1.2'] = 'MHb_C'
multiome_sce$refined_mid_cluster[multiome_sce$refined_mid_cluster == 'MHb.3'] = 'MHb_D'
multiome_sce$refined_mid_cluster[multiome_sce$refined_mid_cluster == 'LHb.2.7'] = 'LHb_A'
multiome_sce$refined_mid_cluster[multiome_sce$refined_mid_cluster == 'LHb.1.3.4'] = 'LHb_B'
multiome_sce$refined_mid_cluster[multiome_sce$refined_mid_cluster == 'LHb.4'] = 'LHb_C'



#Now add the inhibitory annotations

# Add the putative Inhibitory neuron annotations to the midSeurat object
inhib_barcodes_1 = rownames(inhib_meta)[inhib_meta$refined_mid_cluster == 'Putative_Inhib_LHb_4.1']
inhib_barcodes_2 = rownames(inhib_meta)[inhib_meta$refined_mid_cluster == 'Putative_Inhib_LHb_4.2']

midSeurat$refined_mid_cluster[rownames(midSeurat[[]]) %in% inhib_barcodes_1] = 'GABA_LHb_C.1'
midSeurat$refined_mid_cluster[rownames(midSeurat[[]]) %in% inhib_barcodes_2] = 'GABA_LHb_C.2'

#table(midSeurat$refined_mid_cluster, midSeurat$mid_cluster)

#And the finer resolution annotations
midSeurat$refined_cluster_ann[rownames(midSeurat[[]]) %in% inhib_barcodes_1] = 'GABA_LHb_C.1'
midSeurat$refined_cluster_ann[rownames(midSeurat[[]]) %in% inhib_barcodes_2] = 'GABA_LHb_C.2'

#table(midSeurat$refined_cluster_ann, midSeurat$cluster_ann)

#And then with the SingleCellExperiment object too
multiome_sce$refined_mid_cluster[rownames(colData(multiome_sce)) %in% inhib_barcodes_1] = 'GABA_LHb_C.1'
multiome_sce$refined_mid_cluster[rownames(colData(multiome_sce)) %in% inhib_barcodes_2] = 'GABA_LHb_C.2'

#table(multiome_sce$refined_mid_cluster,multiome_sce$mid_cluster)

#Finer resolution annotations
multiome_sce$refined_cluster_ann[rownames(colData(multiome_sce)) %in% inhib_barcodes_1] = 'GABA_LHb_C.1'
multiome_sce$refined_cluster_ann[rownames(colData(multiome_sce)) %in% inhib_barcodes_2] = 'GABA_LHb_C.2'

#table(multiome_sce$refined_cluster_ann, multiome_sce$cluster_ann)





#Double check
table(multiome_sce$refined_mid_cluster, multiome_sce$refined_cluster_ann)



#Some summary plots
#UMAP with the new annotations
#And the donor proportion barplot per cluster


#colors
source(here('code','05_03_annotation_adjustments','celltype_colors.R'))

# ## assign color gradients to mid resolution clusters based on Broad cell-types

# # extract LHb and MHb clusters
# cluster_levels <- c(levels(midSeurat), 'Ependymal')
# cluster_levels
# LHb_clusters <- grep("LHb", cluster_levels, value = TRUE)
# MHb_clusters <- grep("MHb", cluster_levels, value = TRUE)

# # Create tonal gradients for LHb and MHb
# LHb_colors <- colorspace::sequential_hcl(length(LHb_clusters), h = 210, c = 80, l = c(30, 80))
# MHb_colors <- colorspace::sequential_hcl(length(MHb_clusters), h = 320, c = 80, l = c(30, 80))

# # Build full cluster color map
# my_colors_mid <- setNames(rep("#bdbdbd", length(cluster_levels)), cluster_levels)
# my_colors_mid[LHb_clusters] <- LHb_colors
# my_colors_mid[MHb_clusters] <- MHb_colors

# # assign base color for other types from your existing palette
# for (category in c("Oligo", "Astrocyte", "OPC", "Microglia", "Endo", "Inhib.Thal", "Excit.Thal", "Thal", 'Ependymal')) {
#     matched <- grep(category, cluster_levels, value = TRUE)
#     my_colors_mid[matched] <- my_colors[[gsub("\\.", "_", category)]]
# }

# my_colors_mid["Inhib_LHb_4.1"] <- "#8B0000"  # Dark red
# my_colors_mid["Inhib_LHb_4.2"] <- "#DC143C"  # Crimson red
# #Adjust color for MHb3, too light
# my_colors_mid["MHb.3"] <- "#56204eff" 

# my_colors_mid["Inhib.Thal"] <- "#9a9fe7"

plt1 <- DimPlot(midSeurat, 
                reduction = "wnn.umap",
                group.by = "refined_mid_cluster",
                cols = my_colors_mid) + 
    labs(title = "WNN cell types (refined Mid-resolution)")

plt1

plt2 <- DimPlot(midSeurat, 
                reduction = "umap.integrated",
                group.by = "refined_mid_cluster",
                cols = my_colors_mid) + 
    labs(title = "RNA UMAP (refined Mid-resolution)")

plt2

plt3 <- DimPlot(midSeurat, 
                reduction = "umap.lsi.integrated",
                group.by = "refined_mid_cluster",
                cols = my_colors_mid) + 
    labs(title = "ATAC UMAP (refined Mid-resolution)")

plt3


ggsave(here(plot_path, "WNN_umap_refined_mid_cluster.pdf"), plt1, width = 7, height = 7, device = 'pdf')
ggsave(here(plot_path, "RNA_umap_refined_mid_cluster.pdf"), plt2, width = 7, height = 7, device = 'pdf')
ggsave(here(plot_path, "ATAC_umap_refined_mid_cluster.pdf"), plt3, width = 7, height = 7, device = 'pdf')



#And now donor proportion plots
donor_order <- multiome_sce@colData |>
  as.data.frame() |>
  group_by(orig.ident) |>
  summarise(cell_counts = n()) |> 
  arrange(desc(cell_counts)) |>
  pull(orig.ident)

cellNum_plot <- multiome_sce@colData |>
  as.data.frame() |>
  mutate(orig.ident = factor(orig.ident, levels = donor_order)) |>
  group_by(orig.ident, refined_mid_cluster) |>
  summarise(cell_count = n(), .groups = "drop") |>
  ggplot(aes(x = orig.ident, y = cell_count, fill = refined_mid_cluster)) +
  geom_col() +
  scale_fill_manual(values = my_colors_mid) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Donor", y = "Cell Count", fill = "Refined Mid Cluster") +
  ggtitle("Multiome clusters per donor")

cellNum_plot

PropcellNum_plot <- multiome_sce@colData |>
  as.data.frame() |>
  mutate(orig.ident = factor(orig.ident, levels = donor_order)) |>
  group_by(orig.ident, refined_mid_cluster) |>
  summarise(cell_count = n(), .groups = "drop") |>
  ggplot(aes(x = orig.ident, y = cell_count, fill = refined_mid_cluster)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = my_colors_mid) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Donor", y = "Proportion", fill = "Refined Mid Cluster") +
  ggtitle("Multiome cluster proportions per donor")

PropcellNum_plot



donor_colors = MetBrewer::met.brewer("Hokusai1", length(donor_order))


cluster_order <- multiome_sce@colData |>
  as.data.frame() |>
  group_by(refined_mid_cluster) |>
  summarise(n = n()) |> 
  arrange(desc(n)) |>
  pull(refined_mid_cluster)

donorNum_plot <- multiome_sce@colData |>
  as.data.frame() |>
  mutate(refined_mid_cluster = factor(refined_mid_cluster, levels = cluster_order)) |>
  group_by(refined_mid_cluster, orig.ident) |>
  summarise(cell_count = n(), .groups = "drop") |>
  ggplot(aes(x = refined_mid_cluster, y = cell_count, fill = orig.ident)) +
  geom_col() +
  scale_fill_manual(values = donor_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Refined Mid Cluster", y = "Cell Count", fill = "Donor") +
  ggtitle("Multiome donors per cluster")

donorNum_plot

PropdonorNum_plot <- multiome_sce@colData |>
  as.data.frame() |>
  mutate(refined_mid_cluster = factor(refined_mid_cluster, levels = cluster_order)) |>
  group_by(refined_mid_cluster, orig.ident) |>
  summarise(cell_count = n(), .groups = "drop") |>
  ggplot(aes(x = refined_mid_cluster, y = cell_count, fill = orig.ident)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = donor_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Refined Mid Cluster", y = "Proportion", fill = "Donor") +
  ggtitle("Multiome donors per cluster")

PropdonorNum_plot 


heat_df <- multiome_sce@colData |>
  as.data.frame() |>
  count(refined_mid_cluster, orig.ident, name = "cell_count") |>
  mutate(
    refined_mid_cluster = factor(refined_mid_cluster, levels = cluster_order),
    orig.ident = factor(orig.ident, levels = donor_order)
  ) |>
  tidyr::complete(refined_mid_cluster, orig.ident, fill = list(cell_count = 0)) |>
  mutate(
    # keep your existing rule: <10 shown as grey
    cell_count_plot = if_else(cell_count < 10 & cell_count > 0, NA_real_, as.numeric(cell_count))
  )

special_df <- heat_df |>
  mutate(
    special_class = case_when(
      cell_count == 0 ~ "No cells",
      cell_count < 10 ~ "< 10 cells",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(special_class))


count_heatmap <- ggplot() +
  # special tiles (black + grey) with legend
  geom_tile(
    data = special_df,
    aes(x = orig.ident, y = refined_mid_cluster, fill = special_class),
    color = "black", linewidth = 0.2
  ) +
  scale_fill_manual(
    name = "Can't pseudobulk",
    values = c("No cells" = "black", "< 10 cells" = "grey70"),
    breaks = c("No cells", "< 10 cells")
  ) +
  ggnewscale::new_scale_fill() +
  # regular tiles (>=10) with continuous legend
  geom_tile(
    data = filter(heat_df, cell_count >= 10),
    aes(x = orig.ident, y = refined_mid_cluster, fill = cell_count),
    color = "black", linewidth = 0.2
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    name = "Cell count",
    trans = "sqrt"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "orig.ident",
    y = "refined_mid_cluster",
    title = "Cell counts per cluster × donor"
  )

count_heatmap


cellNum_plot
PropcellNum_plot

donorNum_plot
PropdonorNum_plot 

count_heatmap


ggsave(here(plot_path, "cluster_cell_num_per_donor_barplot.pdf"), cellNum_plot, width = 7, height = 7, device = 'pdf')
ggsave(here(plot_path, "cluster_cell_num_per_donor_proportionbarplot.pdf"), PropcellNum_plot, width = 7, height = 7, device = 'pdf')

ggsave(here(plot_path, "donor_num_per_cluster_barplot.pdf"), donorNum_plot, width = 7, height = 7, device = 'pdf')
ggsave(here(plot_path, "donor_num_per_cluster_proportionbarplot.pdf"), PropdonorNum_plot , width = 7, height = 7, device = 'pdf')

ggsave(here(plot_path, "cell_count_donor_cluster_heatmap.pdf"), count_heatmap, width = 7, height = 9, device = 'pdf')


#########
# Add doublet proportional plots for the redone mid clusters and fine as well
#########

#Fine resolution clusters
cluster_order <- midSeurat@meta.data |>
  dplyr::count(refined_cluster_ann, scDblFinder.class) |>
  dplyr::group_by(refined_cluster_ann) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::filter(scDblFinder.class == "doublet") |>
  dplyr::arrange(prop) |>
  dplyr::pull(refined_cluster_ann)

fineCluster_doublet_class_p = midSeurat@meta.data |>
  dplyr::mutate(refined_cluster_ann = factor(refined_cluster_ann, levels = cluster_order)) |>
  ggplot(aes(x = refined_cluster_ann, fill = scDblFinder.class)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Fine clusters", y = "Proportion", fill = "scDblFinder") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

fineCluster_doublet_score_p = midSeurat@meta.data |>
  dplyr::mutate(refined_cluster_ann = factor(refined_cluster_ann, levels = cluster_order)) |>
  ggplot(aes(x = refined_cluster_ann, y = scDblFinder.score)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Fine clusters") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

fineCluster_doublet_class_p 
fineCluster_doublet_score_p

plt5 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "umap.integrated",
                group.by = "refined_cluster_ann",
                label.size = 3) + 
    NoLegend() +
    labs(title = "RNA UMAP (fine-resolution)")

plt5

plt6 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "umap.integrated",
                group.by = "scDblFinder.class",
                label.size = 3) + 
    NoLegend() +
    labs(title = "RNA UMAP (doublet class)")

plt6


ggsave(here(plot_path, "fineCluster_doublet_class_proportionalbarplot.pdf"), fineCluster_doublet_class_p , width = 7, height = 7, device = 'pdf')
ggsave(here(plot_path, "fineCluster_doublet_score_boxplot.pdf"), fineCluster_doublet_score_p, width = 7, height = 7, device = 'pdf')

ggsave(here(plot_path, "rna_umap_fineCluster_annots.pdf"), plt5, width = 7, height = 7, device = 'pdf')
ggsave(here(plot_path, "rna_umap_fineCluster_doublet_annots.pdf"), plt6, width = 7, height = 7, device = 'pdf')



#And save the singlecellexperiment and seurat objects

qs_save(multiome_sce, paste0(new_data_path, '/refined_annotation_multiomeHab_SCE.qs2'))
qs_save(midSeurat, paste0(new_data_path, '/refined_annotation_multiomeHab_Seurat.qs2'))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
