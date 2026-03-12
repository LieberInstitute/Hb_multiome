
#Sub clustering the LHb.4 population
#The co-GABA-GLut population is actually very small, and this is the largest cluster we have for the habenula.
#See if we can resolve some higher resolution clusters.



library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
library(MetaNeighbor)
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


multiome_path = here('processed-data', '08_spatial_registration_vs_multiome_snRNA-seq','mid')

#Multiome human data
multiome_sce = readRDS(paste0(multiome_path, '/seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v5.rds'))
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))


#Set up color scale
color_palette_1 = c('Astrocyte' = 'grey','Endo' = 'firebrick','Microglia' = 'darkred',
'Oligo' = 'darkgoldenrod','OPC' = 'cornsilk3')

color_palette_2 = c('Inhib.Thal' = 'dodgerblue','Excit.Thal' = 'indianred2','Thal' = 'lightsteelblue')

color_palette_3 = MetBrewer::met.brewer("Redon", n = 10)
names(color_palette_3) = c('LHb.1','LHb.1.3','LHb.1.3.4','LHb.2.7','LHb.4','LHb.7','MHb.1','MHb.1.2','MHb.2','MHb.3')

color_palette = c(color_palette_1, color_palette_2, color_palette_3)



#Here, we're going to start with the default scran and scatter, runPCA, runUMAP functions
#The logcounts assay is the default used for PCA

#Just the medial habenula clusters
Lhab_sce = multiome_sce[, multiome_sce$merged_cluster == 'LHb']
#Lhab_sce = multiome_sce[, multiome_sce$mid_cluster %in% c('LHb.4', 'LHb.7')]


rm(multiome_sce)
gc() 

#This is uncorrected data, across all the donors
# 1. Identify highly variable genes
dec <- modelGeneVar(Lhab_sce )
hvg <- getTopHVGs(dec, n = 2000)

# 2. Run PCA on highly variable genes
Lhab_sce <- runPCA(Lhab_sce, subset_row = hvg)

# 3. Run UMAP on PCA space
Lhab_sce <- runUMAP(Lhab_sce, dimred = "PCA", n_dimred = 30)

# 4. Visualize
plotPCA(Lhab_sce, colour_by = "mid_cluster")
plotUMAP(Lhab_sce, colour_by = "mid_cluster")
plotUMAP(Lhab_sce, colour_by = "orig.ident")
plotUMAP(Lhab_sce, colour_by = "cluster_ann")

#All the cell names are unique, so I can break apart by donor and then put it back together
table(duplicated(rownames(colData(Lhab_sce))))

# Get unique donor IDs
donors <- unique(Lhab_sce$orig.ident)

# Create a named list of SCE objects, one per donor
Lhab_sce_list <- lapply(donors, function(donor) {
  Lhab_sce[, Lhab_sce$orig.ident == donor]
})
names(Lhab_sce_list) <- donors

#Get the markers per donor
donor_markers_list = lapply(donors, function(donor) {
  MetaMarkers::compute_markers(assay(Lhab_sce_list[[donor]], 'cpm'), Lhab_sce_list[[donor]]$cluster_ann)
})
names(donor_markers_list) = donors


#Make metamarkers
cross_donor_hab_markers = make_meta_markers(donor_markers_list, detailed_stats = TRUE)

cross_donor_hab_markers %>% filter(rank <= 20) %>% View()

mc1_pareto = plot_pareto_markers(cross_donor_hab_markers, "C.40.LHb.4", min_recurrence = 0) + ggtitle('C.40.LHb.4')
mc1_pareto



cellNum_plot = Lhab_sce@colData %>%
  as.data.frame() %>%
  group_by(orig.ident, cluster_ann) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(orig.ident = forcats::fct_reorder(orig.ident, total, .desc = TRUE)) %>%
  ggplot(aes(x = orig.ident, y = n, fill = cluster_ann)) +
  geom_col() +
  #scale_fill_manual(values = color_palette) +
  #scale_alpha_manual(values = c("≥10 MHb.3" = 1, "<10 MHb.3" = 0.5), name = "MHb.3 cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Donor", y = "Cell Count", fill = "High Cluster") + ggtitle('Lateral habenula clusters per donor')

cellNum_plot 

#ggsave(filename = 'lateral_hab_cluster_donor_numbers.pdf', path = plot_path, 
#plot = cellNum_plot, device = 'pdf', width = 8, height = 6)


#Donor S04_Hb_r is pretty much entirely LHb.4, throws off the cross-donor comparisons
#Also exclude the 03, 05, 09, unbalanced cell-type proportions compared to the others
exclude_donors = c("S04_Hb_r", "S03_Hb_r", "S05_Hb_r", "S09_Hb_r")




#Convert first to seurats
multiome_seurat_list <- lapply(Lhab_sce_list, function(sce) {
  as.Seurat(sce, counts = "counts", data = "cpm")
})

#Clustering for each individual donor
for (i in seq_along(multiome_seurat_list)) {
  # Find variable features
  multiome_seurat_list[[i]] <- FindVariableFeatures(multiome_seurat_list[[i]], selection.method = "vst", nfeatures = 2000)
  
  # Scale data
  all.genes <- rownames(multiome_seurat_list[[i]])
  multiome_seurat_list[[i]] <- ScaleData(multiome_seurat_list[[i]], features = all.genes)
  
  # PCA
  multiome_seurat_list[[i]] <- RunPCA(multiome_seurat_list[[i]], 
    features = VariableFeatures(object = multiome_seurat_list[[i]]))
  
  # Find neighbors and clusters
  multiome_seurat_list[[i]] <- FindNeighbors(multiome_seurat_list[[i]], dims = 1:20)
  multiome_seurat_list[[i]] <- FindClusters(multiome_seurat_list[[i]], resolution = .8)
  
  # UMAP
  multiome_seurat_list[[i]] <- RunUMAP(multiome_seurat_list[[i]], dims = 1:20)
  
  cat("Finished processing", names(multiome_seurat_list)[i], "\n")
}



#Check out the UMAPs for each donor
for (i in seq_along(multiome_seurat_list)) {
  p <- DimPlot(multiome_seurat_list[[i]], reduction = "umap", group.by = 'seurat_clusters', label = TRUE) + 
    ggtitle(names(multiome_seurat_list)[i])
  print(p)
  p <- DimPlot(multiome_seurat_list[[i]], reduction = "umap",group.by = 'cluster_ann' , label = TRUE) + 
    ggtitle(names(multiome_seurat_list)[i]) # + scale_color_manual(values = color_palette)
  print(p)

}



#Convert to SingleCellExperiment objects for MetaNeighbor
donor_sce_list <- lapply(multiome_seurat_list, as.SingleCellExperiment)
#Change the assays names to counts cpm scaledata
for (i in seq_along(donor_sce_list)) {
  names(assays(donor_sce_list[[i]])) <- c("counts", "cpm", "scaledata" )
}


#And run MetaNeighbor on the default Seurat clusters
#Get single SCE object
all_donor_sce = mergeSCE(donor_sce_list)
View(as.data.frame(colData(all_donor_sce)))

#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(dat = all_donor_sce, min_recurrence = 2, exp_labels = all_donor_sce$study_id)
length(global_hvgs)
keep_global_hvgs = global_hvgs[1:2000]


MN_aurocs = MetaNeighborUS(var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$seurat_clusters,
  fast_version = TRUE)

#Plot allby-all AUROC heatmap
MetaNeighbor::plotHeatmap(MN_aurocs, 
  cex = .5,
  title = "MetaNeighbor AUROCs for HabMulti-ome donors")


#And the best versus next approach

MN_best_aurocs = MetaNeighborUS(var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$seurat_clusters,
  fast_version = TRUE,
one_vs_best = TRUE, symmetric_output = FALSE)

#Plot best_vs_next AUROC heatmap
MetaNeighbor::plotHeatmap(MN_best_aurocs, 
  cex = .5,
  title = "MetaNeighbor best_vs_next AUROCs for HabMulti-ome donors")


#Get the metaclusters from the best vs next results, add those annotations to the full SCE object
mclusters = extractMetaClusters(MN_best_aurocs, threshold = .5)
mclusters


#pdf(file = paste0(plot_path, '/lateral_hab_MetaNeighbor_allbyall_aurocs.pdf'), width = 8, height = 8, useDingbats = FALSE)
#MetaNeighbor::plotHeatmap(MN_aurocs, 
#  cex = .5,
#  title = "MetaNeighbor AUROCs for HabMulti-ome donors")
#dev.off()


#pdf(file = paste0(plot_path, '/lateral_hab_MetaNeighbor_bestbynext_aurocs.pdf'), width = 8, height = 8, useDingbats = FALSE)
#MetaNeighbor::plotHeatmap(MN_best_aurocs, 
#  cex = .5,
#  title = "MetaNeighbor best_vs_next AUROCs for HabMulti-ome donors")
#dev.off()



full_cluster_study_labels = paste(all_donor_sce$study_id, all_donor_sce$seurat_clusters, sep = "|")

# Create a vector of meta_cluster names for each element in mclusters
meta_cluster_names <- rep(names(mclusters), sapply(mclusters, length))
# Flatten mclusters to match the order
flat_mclusters <- unlist(mclusters)
# Create a lookup vector
mclusters_lookup <- setNames(meta_cluster_names, flat_mclusters)
# Map each cell's label to its meta_cluster
all_donor_sce$meta_cluster <- mclusters_lookup[full_cluster_study_labels]


#Confusion matrix between metaclusters and the mid_cluster annotations
all_celltype_conf_mat = as.matrix(table(all_donor_sce$meta_cluster, all_donor_sce$mid_cluster))
all_celltype_sum_vec = colSums(all_celltype_conf_mat)
all_celltype_conf_mat  = sweep(all_celltype_conf_mat , 2, all_celltype_sum_vec, "/")

col_fun = circlize::colorRamp2(c(0, 1), c("white", "red"))
ComplexHeatmap::Heatmap(all_celltype_conf_mat, name = 'Proportion of cells', col = col_fun, column_title = 'Metacluster vs integrated mid-res' ,
cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE )


logcounts(all_donor_sce) = log1p(cpm(all_donor_sce))
dec <- modelGeneVar(all_donor_sce)
hvg <- getTopHVGs(dec, n = 2000)

all_donor_sce <- runPCA(all_donor_sce, subset_row = hvg)
all_donor_sce <- runUMAP(all_donor_sce, dimred = "PCA", n_dimred = 20)

plotUMAP(all_donor_sce, colour_by = "mid_cluster") +
  scale_color_manual(values = color_palette) + ggtitle('Unintegrated mid_res clusters')

plotUMAP(all_donor_sce, colour_by = "meta_cluster") +
  ggtitle('Unintegrated meta clusters')


latHab_donors = names(donor_sce_list)
#Get the markers per donor, for the metaclusters
donor_markers_list = lapply(latHab_donors, function(donor) {
  donor_sce = all_donor_sce[, all_donor_sce$study_id == donor]
  MetaMarkers::compute_markers(assay(donor_sce, 'cpm'), donor_sce$meta_cluster)
})
names(donor_markers_list) = latHab_donors

#Make metamarkers
cross_donor_hab_markers = make_meta_markers(donor_markers_list, detailed_stats = TRUE)

cross_donor_hab_markers %>% filter(rank <= 20) %>% View()

#mc1_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster1", min_recurrence = 0) + ggtitle('MetaCluster 1')
#mc2_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster2", min_recurrence = 0) + ggtitle('MetaCluster 2')
#mc3_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster3", min_recurrence = 0) + ggtitle('MetaCluster 3')
#mc4_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster4", min_recurrence = 0) + ggtitle('MetaCluster 4')
#mc5_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster5", min_recurrence = 0) + ggtitle('MetaCluster 5')
#mc6_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster6", min_recurrence = 0) + ggtitle('MetaCluster 6')
#mc7_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster7", min_recurrence = 0) + ggtitle('MetaCluster 7')
#mc8_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster8", min_recurrence = 0) + ggtitle('MetaCluster 8')
#mc9_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster9", min_recurrence = 0) + ggtitle('MetaCluster 9')

#mc1_pareto
#mc2_pareto
#mc3_pareto
#mc4_pareto
#mc5_pareto
#mc6_pareto
#mc7_pareto
#mc8_pareto
#mc9_pareto



#Try default seurat integration

# Find integration anchors
integration_anchors <- FindIntegrationAnchors(object.list = multiome_seurat_list)

# Integrate the datasets
multiome_seurat_integrated <- IntegrateData(anchorset = integration_anchors)

# Set the integrated assay as active
DefaultAssay(multiome_seurat_integrated) <- "integrated"

# Run standard analysis on integrated data
multiome_seurat_integrated <- ScaleData(multiome_seurat_integrated)
multiome_seurat_integrated <- RunPCA(multiome_seurat_integrated)
multiome_seurat_integrated <- RunUMAP(multiome_seurat_integrated, dims = 1:20)


multiome_seurat_integrated[[]]
full_cluster_study_labels = paste(multiome_seurat_integrated$orig.ident, 
  multiome_seurat_integrated$seurat_clusters, sep = "|")

multiome_seurat_integrated$meta_cluster <- unname(mclusters_lookup[full_cluster_study_labels])


# Visualize
int_author_clust_plot = DimPlot(multiome_seurat_integrated, reduction = "umap", group.by = "cluster_ann", pt.size = 1) + ggtitle('high_cluster')
int_author_midclust_plot = DimPlot(multiome_seurat_integrated, reduction = "umap", group.by = "mid_cluster", pt.size = 1) + 
  scale_color_manual(values = color_palette)
int_meta_clust_plot = DimPlot(multiome_seurat_integrated, reduction = "umap", group.by = "meta_cluster", pt.size = 1)
int_donor_plot = DimPlot(multiome_seurat_integrated, reduction = "umap", group.by = "orig.ident", pt.size = 1)

int_author_clust_plot
int_author_midclust_plot
int_meta_clust_plot
int_donor_plot


multiome_seurat_integrated
DefaultAssay(multiome_seurat_integrated) <- "RNA"
FeaturePlot(multiome_seurat_integrated, features = "GAD2", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red")
FeaturePlot(multiome_seurat_integrated, features = "SLC17A6", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red")



DefaultAssay(multiome_seurat_integrated) <- "RNA"

# Get expression data
umap_data <- as.data.frame(Embeddings(multiome_seurat_integrated, reduction = "umap"))
umap_data$GAD2 <- FetchData(multiome_seurat_integrated, vars = "GAD2", slot = "data")[, 1]
umap_data$SLC17A6 <- FetchData(multiome_seurat_integrated, vars = "SLC17A6", slot = "data")[, 1]
umap_data$SLC32A1 <- FetchData(multiome_seurat_integrated, vars = "SLC32A1", slot = "data")[, 1]

# Create a coexpression category
umap_data$coexpression <- ifelse(umap_data$GAD2 > 0 & umap_data$SLC17A6 > 0, "Both",
                                  ifelse(umap_data$GAD2 > 0, "GAD2 only",
                                         ifelse(umap_data$SLC17A6 > 0, "SLC17A6 only", "Neither")))

ggplot(umap_data, aes(x = umap_1, y = umap_2, color = coexpression)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Both" = "purple", "GAD2 only" = "red", "SLC17A6 only" = "blue", "Neither" = "lightgrey")) +
  theme_bw() +
  labs(title = "GAD2 and SLC17A6 Co-expression")


umap_data$coexpression <- ifelse(umap_data$SLC32A1 > 0 & umap_data$SLC17A6 > 0, "Both",
                                  ifelse(umap_data$SLC32A1 > 0, "SLC32A1 only",
                                         ifelse(umap_data$SLC17A6 > 0, "SLC17A6 only", "Neither")))

ggplot(umap_data, aes(x = umap_1, y = umap_2, color = coexpression)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Both" = "purple", "SLC32A1 only" = "red", "SLC17A6 only" = "blue", "Neither" = "lightgrey")) +
  theme_bw() +
  labs(title = "SLC32A1 and SLC17A6 Co-expression")



