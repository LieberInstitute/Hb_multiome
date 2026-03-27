
#Getting better summary plots for the cross-species metaNeighbor analysis
#Just need the saved auroc matrics and the default MetaNeighbor cluster graph, improve on that one


library(MetaNeighbor)
library(SingleCellExperiment)
library(dplyr)
library(here)


#Path to previous generated cross-species aurocs
auroc_data_path = here('processed-data', '98_external_Hb_comparisons','09_humanMultiome_data')
#Path to plot directory
plot_path = here('plots', '98_external_Hb_comparisons', '10_summary_cross_species_plots')

if (!dir.exists(plot_path)) dir.create(plot_path)


#Aurocs
cross_species_bestVsnext_aurocs = readRDS(file = paste0(auroc_data_path, '/cross_species_MN_BestVsNext_neurons.rds'))

cross_species_bestVsnext_aurocs[1:5, 1:5]


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



