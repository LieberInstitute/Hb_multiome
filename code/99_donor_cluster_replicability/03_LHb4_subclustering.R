
#Previous script we were able to isolate GABAergic clusters from LHb4
#It's still an open question whether LHb4 is habenula versus thalamus
#This cluster spatially is on the border and it has the lowest expression of GPR151 and POU4F1, which are the canonical habenula markers we've used in the past.
#Grab habenula vs thalamus markers excluding the LHb4 and LHb7 clusters andd check out their expression
#



library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
#library(MetaNeighbor)
library(dplyr)
library(ggplot2)
library(scater)
library(scran)
library(here)

here::here()

#Path for new data generated
new_data_path = here('processed-data', '99_donor_cluster_replicability','03_LHb4_subclustering')
#Path to plot directory
plot_path = here('plots', '99_donor_cluster_replicability','03_LHb4_subclustering')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)

#Source the bubble plot functions
source(here('code','98_external_Hb_comparisons', 'bubble_plot_functions.R'))


multiome_path = here('processed-data', '08_spatial_registration_vs_multiome_snRNA-seq','mid')

#Multiome human data
multiome_sce = readRDS(paste0(multiome_path, '/seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v5.rds'))
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))
colnames(colData(multiome_sce))

#This is the donor integrated LHb4 and 7 from the previous script, can use to visualize marker expression
script_02_data_path = here('processed-data', '99_donor_cluster_replicability','02_LHb4_investigation')
multiome_seurat_integrated = readRDS(paste0(script_02_data_path, '/multiome_LHb4_LHb7_integrated_seurat.rds'))
colnames(multiome_seurat_integrated[[]])


#Check any previously published habenula markers, from the Yalcinbas paper
#POU4F1, GPR151, CHRNB4, HTR2C, Habenula markers
#LYPD6B, ADARB2, RORB, Thalamus markers

p_bubble = get_bubble_plot_sce(multiome_sce, 
  top_markers = c('POU4F1', 'GPR151', 'HTR2C', 'GAP43','SNAP25', 
  'LYPD6B', 'ADARB2', 'RORB'),
 sample_name = "Multiome Habenula", group_col = "mid_cluster")


#Low levels of expression of thalamas markers within LHb4 and 7
p_bubble[[1]]
p_bubble[[2]]


#Grab those habenula vs thalamus markers excluding LHB4 and 7
not_lhb47_sce = multiome_sce[,!multiome_sce$mid_cluster %in% c('LHb.2.7')]
#Group the habenula clusters
not_lhb47_sce$large_clusters = not_lhb47_sce$merged_cluster
not_lhb47_sce$large_clusters[not_lhb47_sce$large_clusters %in% c('MHb', 'LHb')] = 'Habenula'
markers_excluded_lhb47 = compute_markers(assay(not_lhb47_sce, 'cpm'), 
not_lhb47_sce$large_clusters)

#Check out top markers
markers_excluded_lhb47 %>% group_by(cell_type) %>% slice_max(auroc, n = 10) %>% View()


#Grab the top ten habenula and thalamus markers and bubble plot them
top_markers_exclude = markers_excluded_lhb47 %>% filter(cell_type %in% c('Habenula','Excit_Thal', 'Inhib_Thal')) %>% 
  group_by(cell_type) %>% slice_max(auroc, n = 20) %>% pull(gene) %>% unique()

p_bubble = get_bubble_plot_sce(multiome_sce, 
  top_markers = top_markers_exclude,
 sample_name = "Multiome Habenula", group_col = "mid_cluster")


#LHb4 makes up most of the lateral hab, definetely has more expression of the LHb markers than thalamus
p_bubble[[1]]
p_bubble[[2]]



#Explore hab and thalamus expression within the integrated LHb4 and 7 seurat object
FeaturePlot(multiome_seurat_integrated, features = "RORB", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'CPM')

DimPlot(multiome_seurat_integrated, group.by = "seurat_clusters", label = TRUE, pt.size = 1)

multiome_seurat_integrated[[]]


#Using an aggregate marker expression approach, quantify how well the thalamus markers compared to the habenula markers perform in predicting the LHb4 cluster

#Here, the aggregate CPMs of a set of marker genes is the scores, the LHb.4 and .7 annotations are the labels

# multiome_sce

# compute_auc = function(scores, label){
#   label = as.logical(label)
#   n1 = as.numeric(sum(label))
#   n2 = as.numeric(sum(!label))
#   R1 = sum(rank(scores)[label])
#   U1 = R1 - n1 * (n1 + 1)/2
#   auc = U1/(n1 * n2)
#   return(auc)
# }

# annot_label = multiome_sce$mid_cluster %in% c('LHb.4', 'LHb.7' )

# predict_lhb47_mat = matrix(NA, nrow = 100, ncol = 3, 
#   dimnames = list(paste0("Top_", 1:100, "_markers"), c('EThal', 'IThal', 'Hab')))

# for(i in 1:100){

#   num_top = i

#   EThal_markers_to_use = markers_excluded_lhb47 %>% filter(cell_type == 'Excit_Thal') %>% 
#     slice_max(auroc, n = num_top) %>% pull(gene)

#   IThal_markers_to_use = markers_excluded_lhb47 %>% filter(cell_type == 'Inhib_Thal') %>% 
#     slice_max(auroc, n = num_top) %>% pull(gene)

#   Hab_markers_to_use = markers_excluded_lhb47 %>% filter(cell_type == 'Habenula') %>% 
#     slice_max(auroc, n = num_top) %>% pull(gene)

#   EThal_agg_expression = colSums(assay(multiome_sce[EThal_markers_to_use, ], 'logcounts'))
#   IThal_agg_expression = colSums(assay(multiome_sce[IThal_markers_to_use, ], 'logcounts'))
#   Hab_agg_expression = colSums(assay(multiome_sce[Hab_markers_to_use, ], 'logcounts'))

#   predict_lhb47_mat[ i, 'EThal'] = compute_auc(EThal_agg_expression, annot_label)
#   predict_lhb47_mat[ i, 'IThal'] = compute_auc(IThal_agg_expression, annot_label)
#   predict_lhb47_mat[ i, 'Hab'] = compute_auc(Hab_agg_expression, annot_label)

# }

# predict_lhb47_mat



# #Plot lines
# plot_df = as.data.frame(predict_lhb47_mat) %>% 
#   tibble::rownames_to_column(var = 'marker_set') %>% 
#   tidyr::pivot_longer(cols = -marker_set, names_to = 'cell_type', values_to = 'auroc') %>%
#   mutate(marker_set = factor(marker_set, levels = paste0("Top_", 1:100, "_markers")))

# ggplot(plot_df, aes(x = marker_set, y = auroc, color = cell_type, group = cell_type)) +
#   geom_line() +
#   theme_bw() +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   labs(x = 'Number of top markers used', y = 'AUROC', color = 'Marker Cell Type') +
#   ggtitle('Predicting LHb4 and LHb7 using aggregate marker expression')



#Also good to visualize in the original multiome UMAP

####################
#Adapting code from code/05_Clustering_ARCr/22_add_mid_level_clustering.R
#https://github.com/LieberInstitute/Hb_multiome/blob/ba13de14c5616b0488eaaab1800199209f89a882/code/05_Clustering_ARCr/22_add_mid_level_clustering.R#L115-L191


## Seurat object with the mid-level cluster annots and dim reductions saved
midSeurat_Dir <- here(
    "processed-data",
    "05_Clustering_ARCr",
    "22_add_mid_level_clustering"
)
Seurat_base_name <- "seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_HD.rds"

mid_file_name <- here(midSeurat_Dir, Seurat_base_name)

midSeurat = readRDS(mid_file_name)

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

Reductions(midSeurat)

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

FeaturePlot(midSeurat, features = "HTR2C", reduction = "wnn.umap", slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'logcounts')



#Data here is log counts
rownames(midSeurat@assays$RNA@layers$data) = rownames(midSeurat@assays$RNA@features)


#Compute markers excluding a cell-type
exclude_celltypes = c('LHb.4')
exclude_celltype_sce = multiome_sce[,!multiome_sce$mid_cluster %in% exclude_celltypes]
#Group the habenula clusters and compute markers
exclude_celltype_sce$large_clusters = exclude_celltype_sce$merged_cluster
exclude_celltype_sce$large_clusters[exclude_celltype_sce$large_clusters %in% c('MHb', 'LHb')] = 'Habenula'
markers_excluded = compute_markers(assay(exclude_celltype_sce, 'cpm'), 
exclude_celltype_sce$large_clusters)


for(i in c(10, 20, 50, 100)){

  num_top = i
  EThal_markers_to_use = markers_excluded %>% filter(cell_type == 'Excit_Thal') %>% 
    slice_max(auroc, n = num_top) %>% pull(gene)

  IThal_markers_to_use = markers_excluded %>% filter(cell_type == 'Inhib_Thal') %>% 
    slice_max(auroc, n = num_top) %>% pull(gene)

  Hab_markers_to_use = markers_excluded %>% filter(cell_type == 'Habenula') %>% 
    slice_max(auroc, n = num_top) %>% pull(gene)

  EThal_agg_expression = colSums(midSeurat@assays$RNA@layers$data[EThal_markers_to_use, ])
  IThal_agg_expression = colSums(midSeurat@assays$RNA@layers$data[IThal_markers_to_use, ])
  Hab_agg_expression = colSums(midSeurat@assays$RNA@layers$data[Hab_markers_to_use, ])

  midSeurat$EThal_agg = EThal_agg_expression / max(EThal_agg_expression)
  midSeurat$IThal_agg = IThal_agg_expression / max(IThal_agg_expression)
  midSeurat$Hab_agg = Hab_agg_expression / max(Hab_agg_expression)

  p_thal = FeaturePlot(midSeurat, features = "EThal_agg", reduction = "wnn.umap") +
    scale_color_gradient(low = "white", high = "red", name = 'aggregate_exp') + 
    ggtitle(paste0("Top ", num_top, " Excit_Thal markers excluding ", exclude_celltypes))

  p_inhib_thal = FeaturePlot(midSeurat, features = "IThal_agg", reduction = "wnn.umap") +
    scale_color_gradient(low = "white", high = "red", name = 'aggregate_exp') + 
    ggtitle(paste0("Top ", num_top, " Inhib_Thal markers excluding ", exclude_celltypes))

  p_hab = FeaturePlot(midSeurat, features = "Hab_agg", reduction = "wnn.umap") +
    scale_color_gradient(low = "white", high = "red", name = 'aggregate_exp') + 
    ggtitle(paste0("Top ", num_top, " Habenula markers excluding ", exclude_celltypes))

  print(p_thal)
  print(p_inhib_thal)
  print(p_hab)
}


#Get the violins of average expression per cluster of the top 50 genes
num_top = 50
EThal_markers_to_use = markers_excluded %>% filter(cell_type == 'Excit_Thal') %>% 
  slice_max(auroc, n = num_top) %>% pull(gene)

IThal_markers_to_use = markers_excluded %>% filter(cell_type == 'Inhib_Thal') %>% 
  slice_max(auroc, n = num_top) %>% pull(gene)

Hab_markers_to_use = markers_excluded %>% filter(cell_type == 'Habenula') %>% 
  slice_max(auroc, n = num_top) %>% pull(gene)




avg_inhib_thal_expr = colMeans(assay(multiome_sce, 'cpm')[IThal_markers_to_use, ])
avg_excite_thal_expr = colMeans(assay(multiome_sce, 'cpm')[EThal_markers_to_use, ])
avg_hab_expr = colMeans(assay(multiome_sce, 'cpm')[Hab_markers_to_use, ])

labels = multiome_sce$mid_cluster

marker_df = data.frame(avg_inhib_thal_expr = avg_inhib_thal_expr, avg_excite_thal_expr = avg_excite_thal_expr, 
  avg_hab_expr = avg_hab_expr, 
  cluster = labels)

#Inhibitory thalamus marker expression
cluster_order <- marker_df |>
  summarize(med = median(avg_inhib_thal_expr, na.rm = TRUE), .by = cluster) |>
  arrange(med) |>
  pull(cluster)

marker_df <- marker_df |>
  mutate(cluster = factor(cluster, levels = cluster_order))

p_avg_inh_thal = ggplot(marker_df, aes(x = cluster, y = avg_inhib_thal_expr, fill = labels)) +
  geom_violin(scale = 'width') +
  theme_bw() +
  scale_fill_manual(values = my_colors_mid) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Cluster") +
  ylab(paste0("Average expression of top ",num_top," inhibitory thalamus markers"))+
  ggtitle(paste0("Inhib. Thalamus marker expression excluding ", exclude_celltypes))


#Excitatory thalamus marker expression
cluster_order <- marker_df |>
  summarize(med = median(avg_excite_thal_expr, na.rm = TRUE), .by = cluster) |>
  arrange(med) |>
  pull(cluster)

marker_df <- marker_df |>
  mutate(cluster = factor(cluster, levels = cluster_order))

p_avg_ect_thal = ggplot(marker_df, aes(x = cluster, y = avg_excite_thal_expr, fill = labels)) +
  geom_violin(scale = 'width') +
  theme_bw() +
  scale_fill_manual(values = my_colors_mid) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Cluster") +
  ylab(paste0("Average expression of top ",num_top," Excitatory thalamus markers"))+
  ggtitle(paste0("Excite. Thalamus marker expression excluding ", exclude_celltypes))

#Habenula marker expression
cluster_order <- marker_df |>
  summarize(med = median(avg_hab_expr, na.rm = TRUE), .by = cluster) |>
  arrange(med) |>
  pull(cluster)

marker_df <- marker_df |>
  mutate(cluster = factor(cluster, levels = cluster_order))

p_avg_hab = ggplot(marker_df, aes(x = cluster, y = avg_hab_expr, fill = labels)) +
  geom_violin(scale = 'width') +
  theme_bw() +
  scale_fill_manual(values = my_colors_mid) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Cluster") +
  ylab(paste0("Average expression of top ",num_top," Habenula markers"))+
  ggtitle(paste0("Habenula marker expression excluding ", exclude_celltypes))

p_avg_hab
p_avg_inh_thal 
p_avg_ect_thal 




#Checking out expression just within the integrated LHb4 clusters
#Check out aggregate expression in the UMAP
num_top = 20
EThal_markers_to_use = markers_excluded %>% filter(cell_type == 'Excit_Thal') %>% 
  slice_max(auroc, n = num_top) %>% pull(gene)

IThal_markers_to_use = markers_excluded %>% filter(cell_type == 'Inhib_Thal') %>% 
  slice_max(auroc, n = num_top) %>% pull(gene)

Hab_markers_to_use = markers_excluded %>% filter(cell_type == 'Habenula') %>% 
  slice_max(auroc, n = num_top) %>% pull(gene)

EThal_agg_expression = colMeans(multiome_seurat_integrated@assays$RNA@data[EThal_markers_to_use, ])
IThal_agg_expression = colMeans(multiome_seurat_integrated@assays$RNA@data[IThal_markers_to_use, ])
Hab_agg_expression = colMeans(multiome_seurat_integrated@assays$RNA@data[Hab_markers_to_use, ])


multiome_seurat_integrated$EThal_agg = EThal_agg_expression
multiome_seurat_integrated$IThal_agg = IThal_agg_expression 
multiome_seurat_integrated$Hab_agg = Hab_agg_expression 

max( EThal_agg_expression, max(IThal_agg_expression , Hab_agg_expression ))

FeaturePlot(multiome_seurat_integrated, features = "EThal_agg", reduction = "umap", pt.size = 1) +
  scale_color_gradient(low = "white", high = "red", name = 'Avg CPM', limits = c(0,4500))

FeaturePlot(multiome_seurat_integrated, features = "IThal_agg", reduction = "umap", pt.size = 1) +
  scale_color_gradient(low = "white", high = "red", name = 'Avg CPM', limits = c(0,4500))

FeaturePlot(multiome_seurat_integrated, features = "Hab_agg", reduction = "umap", pt.size = 1) +
  scale_color_gradient(low = "white", high = "red", name = 'Avg CPM', limits = c(0,4500))

DimPlot(multiome_seurat_integrated, group.by = "seurat_clusters", label = TRUE, pt.size = 1)

FeaturePlot(multiome_seurat_integrated, features = "GPR151", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'CPM')

FeaturePlot(multiome_seurat_integrated, features = "POU4F1", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'CPM')

FeaturePlot(multiome_seurat_integrated, features = "HTR2C", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'CPM')

FeaturePlot(multiome_seurat_integrated, features = "RBFOX1", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'CPM')

FeaturePlot(multiome_seurat_integrated, features = "OTX2", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'CPM')

FeaturePlot(multiome_seurat_integrated, features = "OTX2-AS1", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'CPM')
