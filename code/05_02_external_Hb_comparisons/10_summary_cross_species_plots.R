
#Getting better summary plots for the cross-species metaNeighbor analysis
#Just need the saved auroc matrics and the default MetaNeighbor cluster graph, improve on that one


library(MetaNeighbor)
library(SingleCellExperiment)
library(dplyr)
library(here)


#Path to previous generated cross-species aurocs
auroc_data_path = here('processed-data', '05_02_external_Hb_comparisons','09_humanMultiome_data')
#Path to plot directory
plot_path = here('plots', '05_02_external_Hb_comparisons', '10_summary_cross_species_plots')

if (!dir.exists(plot_path)) dir.create(plot_path)


#Aurocs
cross_species_bestVsnext_aurocs = readRDS(file = paste0(auroc_data_path, '/cross_species_MN_BestVsNext_neurons.rds'))

cross_species_bestVsnext_aurocs[1:5, 1:5]

cross_species_allVsAll_aurocs = readRDS(file = paste0(auroc_data_path, '/cross_species_MN_AllVsAll_neurons.rds'))
cross_species_allVsAll_aurocs[1:5, 1:5]

#merged SCE object
all_donor_sce = readRDS( file = paste0(auroc_data_path, '/cross_species_merged_sce_neurons.rds'))
all_donor_sce

cross_species_bestVsnext_aurocs[1:10, 1:10]
cluster_graph = makeClusterGraph(cross_species_bestVsnext_aurocs , low_threshold = .5)

plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)

#Make the default graph a bit nicer
#Species by color
#Dataset gets its own shape

metadata_df = as.data.frame(colData(all_donor_sce))
metadata_df$species <- gsub("\\|", ".", metadata_df$species)
metadata_df %>% View()


temp_cluster_graph = cluster_graph

species_color_palette = c('Wallace.Mouse' = 'brown','Hashikawa.Mouse' = 'brown',
 'yalcinbas.Human' = 'purple', 'multiome.Human' = 'purple', 'Zebrafish' = 'darkblue')

species_shape_palette = c('Wallace.Mouse' = 'square',
                          'Hashikawa.Mouse' = 'circle',
                          'yalcinbas.Human' = 'square',
                          'multiome.Human' = 'circle',
                          'Zebrafish' = 'circle')

#Subclass colors for the mouse clusters
species_colors = species_color_palette[metadata_df$species]
names(species_colors) = paste(metadata_df$species, metadata_df$final_Annotations, sep = '|')

vertex_study <- factor(igraph::V(temp_cluster_graph)$name)

#Set up the vertex colors
levels(vertex_study) <- species_colors[levels(vertex_study)]

igraph::V(temp_cluster_graph)$color <- as.character(vertex_study)
igraph::V(temp_cluster_graph)$label.color <- "black"


# Create shape vector similar to species_colors
species_shapes = species_shape_palette[metadata_df$species]
names(species_shapes) = paste(metadata_df$species, metadata_df$final_Annotations, sep = '|')

# Set up the vertex shapes
vertex_study_shape <- factor(igraph::V(temp_cluster_graph)$name)
levels(vertex_study_shape) <- species_shapes[levels(vertex_study_shape)]

igraph::V(temp_cluster_graph)$shape <- as.character(vertex_study_shape)


#Just the cluster name as the label, excluding the study
full_vertex_labels = igraph::V(temp_cluster_graph)$name
vertex_labels = sapply(strsplit(igraph::V(temp_cluster_graph)$name, split = '|', fixed = T), '[[', 2)
igraph::V(temp_cluster_graph)$label <- vertex_labels


igraph::E(temp_cluster_graph)$width <- igraph::E(temp_cluster_graph)$weight
igraph::E(temp_cluster_graph)$color <- c("orange","darkgray")[as.numeric(igraph::E(temp_cluster_graph)$weight >= 0.5) + 1]

size_factor=1
plot(temp_cluster_graph,
     vertex.size = 8,
     vertex.label.cex=1.5,
     vertex.label.family="Helvetica",
     vertex.label.font=1,
     vertex.frame.color = 'black',
     edge.width = igraph::E(temp_cluster_graph)$width * 5,
     edge.arrow.size=.1*5,
     edge.arrow.width=0.5*5)
graphics::legend(
    "topleft", 
    legend = names(species_color_palette), 
    pt.bg = species_color_palette,
    col = species_color_palette,
    pt.cex = 2, cex = 0.5*2, pch=c(15, 16, 15, 16, 16)
)

pdf(paste0(plot_path, '/clean_cluster_graph_low_MN_human_mouse_zeb_neurons.pdf'), width = 10, height = 8)
plot(temp_cluster_graph,
     vertex.size = 8,
     vertex.label.cex=1.5,
     vertex.label.family="Helvetica",
     vertex.label.font=1,
     vertex.frame.color = 'black',
     edge.width = igraph::E(temp_cluster_graph)$width * 5,
     edge.arrow.size=.1*5,
     edge.arrow.width=0.5*5)
graphics::legend(
    "topleft", 
    legend = names(species_color_palette), 
    pt.bg = species_color_palette,
    col = species_color_palette,
    pt.cex = 2, cex = 0.5*2, pch=c(15, 16, 15, 16, 16)
)
dev.off()



#And replot the all by all heatmap, adding a species annotation to the heatmap
#Will likely need to grab the MetaNeighbor code and manually adjust


#double check this is the matrix as expected
plotHeatmap(
  cross_species_allVsAll_aurocs,
  cex = .5
)
title("MetaNeighbor Human, Mouse, ZebF Habenula: 2000 HVGs")



#Adapting Original MetaNeighbor code from https://github.com/gillislab/MetaNeighbor/blob/a0aad62868621241caef857350318ff674d94704/R/visualization.R
orderCellTypes <- function(M, na_value = 0) {
    M <- (M + t(M))/2
    M[is.na(M)] <- na_value
    result <- stats::hclust(stats::as.dist(1-M), method = "average")
    return(result)
}

# Extract species from rownames (text before the '|'), keeping only species name
species_labels <- sub("\\|.*", "", rownames(cross_species_allVsAll_aurocs))
species_labels <- gsub("Wallace\\.|Hashikawa\\.|yalcinbas\\.|multiome\\.", "", species_labels)

species_colors <- c(
  "Mouse"     = "#a72a2a",
  "Zebrafish" = "#2a2d7c",
  "Human"     = "#8151a1"
)

species_row_anno <- ComplexHeatmap::HeatmapAnnotation(
  Species = species_labels,
  col = list(Species = species_colors),
  which = "row",
  annotation_name_gp = grid::gpar(fontsize = 8),
  simple_anno_size = grid::unit(4, "mm")
)

annot_heatmap = ComplexHeatmap::Heatmap(cross_species_allVsAll_aurocs,
  name = "AUROC",
  col = rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100)),
  na_col = "gray90",
  cluster_rows = orderCellTypes(cross_species_allVsAll_aurocs),
  cluster_columns = orderCellTypes(cross_species_allVsAll_aurocs),
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 8),
  column_names_gp = grid::gpar(fontsize = 8),
  right_annotation = species_row_anno,
  row_dend_width = grid::unit(15, "mm"),
  column_dend_height = grid::unit(15, "mm")
)

annot_heatmap = ComplexHeatmap::draw(annot_heatmap)


pdf(paste0(plot_path, '/AllvsAll_MN_human_mouse_zeb_neurons.pdf'), width = 10, height = 8)
ComplexHeatmap::draw(annot_heatmap)
dev.off()