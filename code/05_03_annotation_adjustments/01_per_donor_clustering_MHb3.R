
#One of the main conclusions from the cross-species comparisons, is that the current MHb.3 population may be human specific
#It is a small cluster, but almost equally contributed to across donors, and does not have a strong match to mouse or zebrafish clusters.

#This analysis step will do a bit more digging into the MHb.3 cluster.
#Central idea, is to see if this cluster arises obviously when clustering at a per-donor level.

#Try subclustering just the medial habenula clusters

library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
library(MetaNeighbor)
library(dplyr)
library(ggplot2)
library(scater)
library(scran)
library(qs2)
library(here)

here::here()

#Path for new data generated
new_data_path = here('processed-data', '99_donor_cluster_replicability','01_per_donor_clustering_MHb3')
#Path to plot directory
plot_path = here('plots', '99_donor_cluster_replicability','01_per_donor_clustering_MHb3')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '05_5_drop_doublets','01_drop_doublets_and_reDimReduce')

#Multiome human data
multiome_sce = qs_read(paste0(multiome_path, '/reprocessed_doubletRemoved_multiomeHab_SCE.qs2'))
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
Mhab_clusters = c('MHb.1','MHb.1.2','MHb.2','MHb.3')
Mhab_sce = multiome_sce[, multiome_sce$merged_cluster == 'MHb']

#rm(multiome_sce)
#gc() 

#This is uncorrected data, across all the donors
# 1. Identify highly variable genes
dec <- modelGeneVar(Mhab_sce )
hvg <- getTopHVGs(dec, n = 2000)

# 2. Run PCA on highly variable genes
Mhab_sce <- runPCA(Mhab_sce, subset_row = hvg)

# 3. Run UMAP on PCA space
Mhab_sce <- runUMAP(Mhab_sce, dimred = "PCA", n_dimred = 30)

# 4. Visualize
plotPCA(Mhab_sce, colour_by = "mid_cluster")
plotUMAP(Mhab_sce, colour_by = "mid_cluster")
plotUMAP(Mhab_sce, colour_by = "orig.ident")

table(Mhab_sce$orig.ident, Mhab_sce$mid_cluster)


#All the cell names are unique, so I can break apart by donor and then put it back together
table(duplicated(rownames(colData(Mhab_sce))))

# Get unique donor IDs
donors <- unique(Mhab_sce$orig.ident)

# Create a named list of SCE objects, one per donor
Mhab_sce_list <- lapply(donors, function(donor) {
  Mhab_sce[, Mhab_sce$orig.ident == donor]
})
names(Mhab_sce_list) <- donors

#Get the markers per donor
donor_markers_list = lapply(donors, function(donor) {
  MetaMarkers::compute_markers(assay(Mhab_sce_list[[donor]], 'cpm'), Mhab_sce_list[[donor]]$mid_cluster)
})
names(donor_markers_list) = donors


table(Mhab_sce$orig.ident, Mhab_sce$mid_cluster)


#test_donor = "S10_Hb_r"

#Donor dimension reduction
#test_sce = Mhab_sce_list[[test_donor]]
#dec <- modelGeneVar(test_sce)
#hvg <- getTopHVGs(dec, n = 2000)

#test_sce <- runPCA(test_sce, subset_row = hvg, exprs_values = "scaledata")
#test_sce <- runUMAP(test_sce, dimred = "PCA", n_dimred = 20)



#Compute metamarkers excluding a donor
#cross_donor_hab_markers = make_meta_markers(donor_markers_list[names(donor_markers_list) != test_donor], detailed_stats = TRUE)

#Project the metamarkers onto the excluded donor
#top_markers = cross_donor_hab_markers %>% filter(cell_type != 'LHb.7' & rank <= 200)
#ct_scores = score_cells(log1p(cpm(Mhab_sce_list[[test_donor]])), top_markers)
#ct_enrichment = compute_marker_enrichment(ct_scores)
#ct_pred = assign_cells(ct_scores)

#test_sce$cross_donor_pred_200 = ct_pred$predicted

#top_markers = cross_donor_hab_markers %>% filter(cell_type != 'LHb.7' & rank <= 100)
#ct_scores = score_cells(log1p(cpm(Mhab_sce_list[[test_donor]])), top_markers)
#ct_enrichment = compute_marker_enrichment(ct_scores)
#ct_pred = assign_cells(ct_scores)

#test_sce$cross_donor_pred_100 = ct_pred$predicted

#top_markers = cross_donor_hab_markers %>% filter(cell_type != 'LHb.7' & rank <= 50) 
#ct_scores = score_cells(log1p(cpm(Mhab_sce_list[[test_donor]])), top_markers)
#ct_enrichment = compute_marker_enrichment(ct_scores)
#ct_pred = assign_cells(ct_scores)

#test_sce$cross_donor_pred_50 = ct_pred$predicted


#plotUMAP(test_sce, colour_by = "mid_cluster") +
#  scale_color_manual(values = color_palette) + ggtitle('Original cluster annots')

#for (cluster in mid_res_annots) {
#  test_sce$highlight <- ifelse(test_sce$mid_cluster == cluster, cluster, "Other")
#  
#  p <- plotUMAP(test_sce, colour_by = "highlight") +
#    scale_color_manual(values = rlang::list2(!!cluster := color_palette[cluster], Other = "lightgrey")) +
#    ggtitle(cluster)
#  
#  print(p)
#}


#plotUMAP(test_sce, colour_by = "cross_donor_pred_200") +
#  scale_color_manual(values = color_palette) + ggtitle('Predicted cross-donor: top 200 markers')

#plotUMAP(test_sce, colour_by = "cross_donor_pred_100") +
#  scale_color_manual(values = color_palette) + ggtitle('Predicted cross-donor: top 100 markers')

#plotUMAP(test_sce, colour_by = "cross_donor_pred_50") +
#  scale_color_manual(values = color_palette) + ggtitle('Predicted cross-donor: top 50 markers')


#table(test_sce$cross_donor_pred_100)



#Just going to go ahead with the cross-donor replicability clustering, just to compare to the integrated cluster labels

#Convert all the donor SCEs to seurats, do default seurat clustering, convert back to SCEs, and then metaneighbor across donors

#Focus on the donors with at least 10 MHb3 cells, the other donors barely have any medial Hab anyhow
table(Mhab_sce$orig.ident, Mhab_sce$mid_cluster)

#Make a plot showing the cell numbers and annotated cell types per cluster


cellNum_plot = Mhab_sce@colData %>%
  as.data.frame() %>%
  group_by(orig.ident, mid_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(total = sum(n),
         mhb3_count = sum(n[mid_cluster == "MHb.3"], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(orig.ident = forcats::fct_reorder(orig.ident, total, .desc = TRUE),
         highlight = ifelse(mhb3_count >= 10, "≥10 MHb.3", "<10 MHb.3")) %>%
  ggplot(aes(x = orig.ident, y = n, fill = mid_cluster, alpha = highlight)) +
  geom_col() +
  scale_fill_manual(values = color_palette) +
  scale_alpha_manual(values = c("≥10 MHb.3" = 1, "<10 MHb.3" = 0.5), name = "MHb.3 cells") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Donor", y = "Cell Count", fill = "Mid Cluster") + ggtitle('Medial habenula clusters per donor')

cellNum_plot 
ggsave(filename = 'medial_hab_cluster_donor_numbers.pdf', path = plot_path, 
plot = cellNum_plot, device = 'pdf', width = 8, height = 6)


keep_donors = c("S03_Hb_r", "S07_Hb_r", "S08_Hb_r", "S10_Hb_r", "S11_Hb_r", 'S12_Hb_r')


#Convert first to seurats
multiome_seurat_list <- lapply(Mhab_sce_list[keep_donors], function(sce) {
  altExps(sce) <- NULL
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

#For trying different clustering resolutions
#for (i in seq_along(multiome_seurat_list)) {
#  multiome_seurat_list[[i]] <- FindClusters(multiome_seurat_list[[i]], resolution = .8)
#}



#Check out the UMAPs for each donor
for (i in seq_along(multiome_seurat_list)) {
  p <- DimPlot(multiome_seurat_list[[i]], reduction = "umap", group.by = 'seurat_clusters', label = TRUE, pt.size = 2) + 
    ggtitle(names(multiome_seurat_list)[i])
  print(p)
  p <- DimPlot(multiome_seurat_list[[i]], reduction = "umap",group.by = 'mid_cluster' , label = TRUE, pt.size = 2) + 
    ggtitle(names(multiome_seurat_list)[i]) + scale_color_manual(values = color_palette)
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

pdf(file = paste0(plot_path, '/medial_hab_MetaNeighbor_allbyall_aurocs.pdf'), width = 8, height = 8, useDingbats = FALSE)
MetaNeighbor::plotHeatmap(MN_aurocs, 
  cex = .5,
  title = "MetaNeighbor AUROCs for HabMulti-ome donors")
dev.off()


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

pdf(file = paste0(plot_path, '/medial_hab_MetaNeighbor_bestbynext_aurocs.pdf'), width = 8, height = 8, useDingbats = FALSE)
MetaNeighbor::plotHeatmap(MN_best_aurocs, 
  cex = .5,
  title = "MetaNeighbor best_vs_next AUROCs for HabMulti-ome donors")
dev.off()

#Get the metaclusters from the best vs next results, add those annotations to the full SCE object
mclusters = extractMetaClusters(MN_best_aurocs, threshold = .6)
mclusters

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

pdf(file = paste0(plot_path, '/confusMat_metacluster_author_annots.pdf'), width = 8, height = 8, useDingbats = FALSE)
ComplexHeatmap::Heatmap(all_celltype_conf_mat, name = 'Proportion of cells', col = col_fun, column_title = 'Metacluster vs integrated mid-res' ,
cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE )
dev.off()

colnames(colData(all_donor_sce))
#Confusion matrix with the higher resolution existing clusters
all_celltype_conf_mat = as.matrix(table(all_donor_sce$meta_cluster, all_donor_sce$cluster_ann))
all_celltype_sum_vec = colSums(all_celltype_conf_mat)
all_celltype_conf_mat  = sweep(all_celltype_conf_mat , 2, all_celltype_sum_vec, "/")

#Filter for just the MHb clusters
all_celltype_conf_mat = all_celltype_conf_mat[, grepl('MHb', colnames(all_celltype_conf_mat ))]

col_fun = circlize::colorRamp2(c(0, 1), c("white", "red"))
ComplexHeatmap::Heatmap(all_celltype_conf_mat, name = 'Proportion of cells', col = col_fun, column_title = 'Metacluster vs integrated high-res' ,
cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE )

pdf(file = paste0(plot_path, '/confusMat_metacluster_highRes_author_annots.pdf'), width = 8, height = 8, useDingbats = FALSE)
ComplexHeatmap::Heatmap(all_celltype_conf_mat, name = 'Proportion of cells', col = col_fun, column_title = 'Metacluster vs integrated high-res' ,
cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE )
dev.off()



logcounts(all_donor_sce) = log1p(cpm(all_donor_sce))
dec <- modelGeneVar(all_donor_sce)
hvg <- getTopHVGs(dec, n = 2000)

all_donor_sce <- runPCA(all_donor_sce, subset_row = hvg)
all_donor_sce <- runUMAP(all_donor_sce, dimred = "PCA", n_dimred = 20)

plotUMAP(all_donor_sce, colour_by = "mid_cluster") +
  scale_color_manual(values = color_palette) + ggtitle('Unintegrated mid_res clusters')

plotUMAP(all_donor_sce, colour_by = "meta_cluster") +
  ggtitle('Unintegrated meta clusters')



medHab_donors = names(donor_sce_list)
#Get the markers per donor, for the metaclusters
donor_markers_list = lapply(medHab_donors, function(donor) {
  donor_sce = all_donor_sce[, all_donor_sce$study_id == donor]
  MetaMarkers::compute_markers(assay(donor_sce, 'cpm'), donor_sce$meta_cluster)
})
names(donor_markers_list) = medHab_donors

#Make metamarkers
cross_donor_hab_markers = make_meta_markers(donor_markers_list, detailed_stats = TRUE)

cross_donor_hab_markers %>% filter(rank <= 20) %>% View()

mc1_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster1", min_recurrence = 0) + ggtitle('MetaCluster 1')
mc2_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster2", min_recurrence = 0) + ggtitle('MetaCluster 2')
mc3_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster3", min_recurrence = 0) + ggtitle('MetaCluster 3')
mc4_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster4", min_recurrence = 0) + ggtitle('MetaCluster 4')
mc5_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster5", min_recurrence = 0) + ggtitle('MetaCluster 5')
mc6_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster6", min_recurrence = 0) + ggtitle('MetaCluster 6')
mc7_pareto = plot_pareto_markers(cross_donor_hab_markers, "meta_cluster7", min_recurrence = 0) + ggtitle('MetaCluster 7')


mc1_pareto
mc2_pareto
mc3_pareto
mc4_pareto
mc5_pareto
mc6_pareto
mc7_pareto

ggsave(plot = mc1_pareto, filename = 'meta_cluster1_pareto_plot.pdf', 
path = plot_path, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)
ggsave(plot = mc2_pareto, filename = 'meta_cluster2_pareto_plot.pdf', 
path = plot_path, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)
ggsave(plot = mc3_pareto, filename = 'meta_cluster3_pareto_plot.pdf', 
path = plot_path, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)
ggsave(plot = mc4_pareto, filename = 'meta_cluster4_pareto_plot.pdf', 
path = plot_path, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)
ggsave(plot = mc5_pareto, filename = 'meta_cluster5_pareto_plot.pdf', 
path = plot_path, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)
ggsave(plot = mc6_pareto, filename = 'meta_cluster6_pareto_plot.pdf', 
path = plot_path, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)
ggsave(plot = mc7_pareto, filename = 'meta_cluster7_pareto_plot.pdf', 
path = plot_path, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)



#From all of the above, annotating the metaclustes
#meta_cluster1 = MHb.1
#meta_cluster2 = MHb.2
#meta_cluster3 = MHb.1.2
#meta_cluster4 = MHb.3
#meta_cluster5 = MHb.1.2
#meta_cluster6 = MHb.1
#meta_cluster7 = oligos



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
int_author_clust_plot = DimPlot(multiome_seurat_integrated, reduction = "umap", group.by = "mid_cluster", pt.size = 1) +
  scale_color_manual(values = color_palette)

int_meta_clust_plot = DimPlot(multiome_seurat_integrated, reduction = "umap", group.by = "meta_cluster", pt.size = 1)
int_donor_plot = DimPlot(multiome_seurat_integrated, reduction = "umap", group.by = "orig.ident", pt.size = 1)

int_author_clust_plot
int_meta_clust_plot
int_donor_plot

ggsave(filename = 'medial_hab_integrated_umap_author_clusters.pdf', path = plot_path, 
plot = int_author_clust_plot, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)

ggsave(filename = 'medial_hab_integrated_umap_meta_clusters.pdf', path = plot_path, 
plot = int_meta_clust_plot, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)

ggsave(filename = 'medial_hab_integrated_umap_donors.pdf', path = plot_path, 
plot = int_donor_plot, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)



#Check out the marker gene panels for the metaclusters
#Source the bubble plot functions
#meta_cluster3 actually has the highest expression of the cholinergic markers, so should be grouped with MHb2
source(here('code','98_external_Hb_comparisons', 'bubble_plot_functions.R'))
DefaultAssay(multiome_seurat_integrated) <- "RNA"
p_bubble = get_bubble_plot(multiome_seurat_integrated, 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1','TAC3', 'GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'MBP'),
 sample_name = "Multiome Medial Habenula", group_col = "meta_cluster")

p_bubble[[1]]
p_bubble[[2]]

pdf(file = paste0(plot_path, '/medial_hab_meta_cluster_bubble_plot_mean_counts.pdf'), width = 10, height = 6, useDingbats = FALSE)
p_bubble[[1]]
dev.off()

pdf(file = paste0(plot_path, '/medial_hab_meta_cluster_bubble_plot_zscore_mean_counts.pdf'), width = 10, height = 6, useDingbats = FALSE)
p_bubble[[2]]
dev.off()



p_bubble = get_bubble_plot(multiome_seurat_integrated, 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1','TAC3', 'GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'MBP'),
 sample_name = "Multiome Medial Habenula", group_col = "mid_cluster")

p_bubble[[1]]
p_bubble[[2]]



#Make adjustments
#meta_cluster4 is MHb.3
#meta_cluster3 is MHb.2
#meta_cluster7 is oligodendrocytes

#Keep everything else the same.

colnames(multiome_seurat_integrated[[]])

multiome_seurat_integrated$refined_mid_cluster = multiome_seurat_integrated$mid_cluster
multiome_seurat_integrated$refined_mid_cluster[multiome_seurat_integrated$meta_cluster == "meta_cluster4"] = "MHb.3"
multiome_seurat_integrated$refined_mid_cluster[multiome_seurat_integrated$meta_cluster == "meta_cluster3"] = "MHb.2"
multiome_seurat_integrated$refined_mid_cluster[multiome_seurat_integrated$meta_cluster == "meta_cluster7"] = "Oligo"

int_author_refined_clust_plot = DimPlot(multiome_seurat_integrated, reduction = "umap", group.by = "refined_mid_cluster", pt.size = 1) +
  scale_color_manual(values = color_palette)
int_author_refined_clust_plot


p_bubble = get_bubble_plot(multiome_seurat_integrated, 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1','TAC3', 'GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'MBP'),
 sample_name = "Multiome Medial Habenula", group_col = "refined_mid_cluster")

p_bubble[[1]]
p_bubble[[2]]

p_bubble = get_bubble_plot(multiome_seurat_integrated, 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1','TAC3', 'GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'MBP'),
 sample_name = "Multiome Medial Habenula", group_col = "cluster_ann")

p_bubble[[1]]
p_bubble[[2]]

#Quick check for the expression of GAD1/2 and VGAT, are there any medial GABAergic cells

med_vgat_p = FeaturePlot(multiome_seurat_integrated, features = "SLC32A1", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'CPM')

med_gad1_p = FeaturePlot(multiome_seurat_integrated, features = "GAD1", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'CPM')

med_gad2_p = FeaturePlot(multiome_seurat_integrated, features = "GAD2", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red", name = 'CPM')

med_vgat_p
med_gad1_p
med_gad2_p

ggsave(filename = 'medial_hab_VGAT_umap.pdf', path = plot_path, 
plot = med_vgat_p, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)

ggsave(filename = 'medial_hab_GAD1_umap.pdf', path = plot_path, 
plot = med_gad1_p, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)

ggsave(filename = 'medial_hab_GAD2_umap.pdf', path = plot_path, 
plot = med_gad2_p, device = 'pdf', width = 8, height = 6, useDingbats = FALSE)


#And check our the updated annotations in the full UMAP

midSeurat = qs_read(paste0(multiome_path, '/reprocessed_doubletRemoved_multiomeHab_seurat.qs2'))

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


# Add the adjusted medial annotations
medial_adjusted_barcodes = rownames(multiome_seurat_integrated[[]])
table(medial_adjusted_barcodes %in% rownames(midSeurat[[]]))

midSeurat$refined_mid_cluster = midSeurat$mid_cluster

index = match(medial_adjusted_barcodes, rownames(midSeurat[[]]))
midSeurat$refined_mid_cluster[index] = multiome_seurat_integrated$refined_mid_cluster

table(midSeurat$refined_mid_cluster, midSeurat$mid_cluster)



plt2 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "refined_mid_cluster",
                label.size = 3,
                cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "WNN cell types (Mid-resolution)")


plt2


#And highlight the changes in red
my_colors_mid['Changed'] <- '#f5093cff'

midSeurat$highlight_annot = midSeurat$mid_cluster
midSeurat$highlight_annot[midSeurat$refined_mid_cluster == 'Oligo' & midSeurat$mid_cluster != 'Oligo'] = 'Changed'

plt3 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "highlight_annot",
                label.size = 3, cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "WNN cell types: updated Oligos in red")


plt3


midSeurat$highlight_annot = midSeurat$mid_cluster
midSeurat$highlight_annot[midSeurat$refined_mid_cluster == 'MHb.2' & midSeurat$mid_cluster != 'MHb.2'] = 'Changed'

plt4 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "highlight_annot",
                label.size = 3, cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "WNN cell types: updated MHb.2 in red")


plt4

midSeurat$highlight_annot = midSeurat$mid_cluster
midSeurat$highlight_annot[midSeurat$refined_mid_cluster == 'MHb.3' & midSeurat$mid_cluster != 'MHb.3'] = 'Changed'

plt5 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "highlight_annot",
                label.size = 3, cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "WNN cell types: updated MHb.3 in red")


plt5


midSeurat$highlight_annot = midSeurat$mid_cluster
midSeurat$highlight_annot[midSeurat$cluster_ann == 'C.11.MHb.1.2' ] = 'Changed'

plt6 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "wnn.umap",
                group.by = "highlight_annot",
                label.size = 3, cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "RNA UMAP cell types: mid resolution")


plt6

plt1
plt2
plt3
plt4
plt5

#And now just the RNA umap with the changed annots in red
midSeurat$highlight_annot = midSeurat$mid_cluster
midSeurat$highlight_annot[midSeurat$refined_mid_cluster == 'Oligo' & midSeurat$mid_cluster != 'Oligo' ] = 'Changed'

plt6 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "umap.integrated",
                group.by = "highlight_annot",
                label.size = 3, cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "RNA UMAP cell types: Updated Oligos in red")


plt6

midSeurat$highlight_annot = midSeurat$mid_cluster
midSeurat$highlight_annot[midSeurat$refined_mid_cluster == 'MHb.2' & midSeurat$mid_cluster != 'MHb.2' ] = 'Changed'

plt7 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "umap.integrated",
                group.by = "highlight_annot",
                label.size = 3, cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "RNA UMAP cell types: Updated MHb2 in red")


plt7


midSeurat$highlight_annot = midSeurat$mid_cluster
midSeurat$highlight_annot[midSeurat$refined_mid_cluster == 'MHb.3' & midSeurat$mid_cluster != 'MHb.3' ] = 'Changed'

plt8 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "umap.integrated",
                group.by = "highlight_annot",
                label.size = 3, cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "RNA UMAP cell types: Updated MHb.3 in red")


plt8

#And now just the C.11.MHb.1.2 cluster in red
midSeurat$highlight_annot = midSeurat$mid_cluster
midSeurat$highlight_annot[midSeurat$cluster_ann == 'C.11.MHb.1.2'  ] = 'Changed'

plt9 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "umap.integrated",
                group.by = "highlight_annot",
                label.size = 3, cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "RNA UMAP cell types: Fine C.11.MHb.1.2 in red")


plt9




FeaturePlot(midSeurat, features = "scDblFinder.score", reduction = "umap.integrated")




#Check out doublet scores of the adjusted annotations. 

midSeurat$new_annot_doublet = midSeurat$mid_cluster
midSeurat$new_annot_doublet[midSeurat$refined_mid_cluster == 'MHb.3' & midSeurat$mid_cluster != 'MHb.3'] = 'New MHb.3'
midSeurat$new_annot_doublet[midSeurat$refined_mid_cluster == 'MHb.2' & midSeurat$mid_cluster != 'MHb.2'] = 'New MHb.2'
midSeurat$new_annot_doublet[midSeurat$refined_mid_cluster == 'Oligo' & midSeurat$mid_cluster != 'Oligo'] = 'New Oligo'

cluster_order <- midSeurat@meta.data |>
  dplyr::count(new_annot_doublet, scDblFinder.class) |>
  dplyr::group_by(new_annot_doublet) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::filter(scDblFinder.class == "doublet") |>
  dplyr::arrange(prop) |>
  dplyr::pull(new_annot_doublet)

midSeurat@meta.data |>
  dplyr::mutate(new_annot_doublet = factor(new_annot_doublet, levels = cluster_order)) |>
  ggplot(aes(x = new_annot_doublet, fill = scDblFinder.class)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Fine clusters", y = "Proportion", fill = "scDblFinder") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

midSeurat@meta.data |>
  dplyr::mutate(new_annot_doublet = factor(new_annot_doublet, levels = cluster_order)) |>
  ggplot(aes(x = new_annot_doublet, y = scDblFinder.score)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Fine clusters") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



#This is the donor integrated LHb4 and 7 from the previous script, can use to visualize marker expression
script_02_data_path = here('processed-data', '99_donor_cluster_replicability','02_LHb4_investigation')
latHb_seurat = readRDS(paste0(script_02_data_path, '/multiome_LHb4_LHb7_integrated_seurat.rds'))
colnames(latHb_seurat[[]])

 
inhib_meta = latHb_seurat[[]]
# Add the putative Inhibitory neuron annotations to the midSeurat object
inhib_barcodes_1 = rownames(inhib_meta)[inhib_meta$refined_mid_cluster == 'Putative_Inhib_LHb_4.1']
inhib_barcodes_2 = rownames(inhib_meta)[inhib_meta$refined_mid_cluster == 'Putative_Inhib_LHb_4.2']

#Exclude them from the rest of the data set

midSeurat$refined_mid_cluster_2 = midSeurat$mid_cluster
midSeurat$refined_mid_cluster_2[rownames(midSeurat[[]]) %in% inhib_barcodes_1] = 'Putative_Inhib_LHb_4.1'
midSeurat$refined_mid_cluster_2[rownames(midSeurat[[]]) %in% inhib_barcodes_2] = 'Putative_Inhib_LHb_4.2'

my_colors_mid["Putative_Inhib_LHb_4.1"] <- "#8B0000"  # Dark red
my_colors_mid["Putative_Inhib_LHb_4.2"] <- "#DC143C"  # Crimson red

plt10 <- DimPlot(midSeurat, 
                label = TRUE, 
                reduction = "umap.integrated",
                group.by = "refined_mid_cluster_2",
                label.size = 3, cols = my_colors_mid) + 
    NoLegend() +
    labs(title = "RNA UMAP cell types")


plt10

