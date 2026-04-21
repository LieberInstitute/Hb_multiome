#Pulling together plots for figures


library(SingleCellExperiment)
library(dplyr)
library(tidyr)
library(ggplot2)
library(qs2)
library(scater)
library(here)

here::here()


#Path to save any generated data
#new_data_path = here('processed-data', '05_03_annotation_adjustments', '09_figure_plots')
#Path to plot directory
plot_path = here('plots','05_03_annotation_adjustments', '09_figure_plots')

#if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))
multiome_sce

SingleCellExperiment::altExp(multiome_sce, "ATAC") <- NULL


source(here('code','05_02_external_Hb_comparisons', 'bubble_plot_functions.R'))


#Colors
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
    Thal = "#4d55b7",
    Ependymal = "#f5a105ff"
)

## assign color gradients to mid resolution clusters based on Broad cell-types

# extract LHb and MHb clusters
cluster_levels <- c(unique(multiome_sce$mid_cluster))
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

my_colors_mid["Inhib_LHb_4.1"] <- "#8B0000"  # Dark red
my_colors_mid["Inhib_LHb_4.2"] <- "#DC143C"  # Crimson red
#Adjust color for MHb3, too light
my_colors_mid["MHb.3"] <- "#56204eff" 
my_colors_mid["Inhib.Thal"] <- "#9a9fe7"
my_colors_mid["Ependymal"] <- "#f5a105ff"

names(my_colors_mid)


#First are the fine to mid cluster adjustments that should have been done before any cross-species
#This is the C.21.Astrocyte being called ependymal and the C.11.MHb.1.2 being called MHb2 based on cholinergic markers

#Show those targeted bubble plots first


get_bubble_plot_sce(multiome_sce,
  top_markers = c('CHAT','SLC18A3','SLC5A7'),
  sample_name = "Multiome Habenula", group_col = "cluster_ann"
)


cell_filt = multiome_sce$mid_cluster %in% c('MHb.3','MHb.2','MHb.1','MHb.1.2')


p_bubble = get_bubble_plot_sce(multiome_sce[ , cell_filt], 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3'),
 sample_name = "Multiome Habenula", group_col = "cluster_ann")

p_bubble[[1]]
p_bubble[[2]]






plt1 <- plotReducedDim(multiome_sce, 
                dimred = "wnn.umap",
                colour_by = "mid_cluster") +
  scale_color_manual(values = my_colors_mid) 

plt1

plt2 <- plotReducedDim(multiome_sce, 
                dimred = "umap.integrated",
                colour_by = "mid_cluster") +
  scale_color_manual(values = my_colors_mid) 

plt2
reducedDimNames(multiome_sce)

table(multiome_sce$cluster_ann)


