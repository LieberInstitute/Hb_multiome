
#One of the main conclusions from the cross-species comparisons, is that the current MHb.3 population may be human specific
#It is a small cluster, but almost equally contributed to across donors, and does not have a strong match to mouse or zebrafish clusters.

#This analysis step will do a bit more digging into the MHb.3 cluster.
#Central idea, is to see if this cluster arises obviously when clustering at a per-donor level.

#Try subclustering just neurons

#Also, check the co-GABA-Glut co-expression in actual cells in LHb.4, but across donors too. 


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
new_data_path = here('processed-data', '99_donor_cluster_replicability','01_per_donor_clustering_MHb3')
#Path to plot directory
plot_path = here('plots', '99_donor_cluster_replicability','01_per_donor_clustering_MHb3')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '08_spatial_registration_vs_multiome_snRNA-seq','mid')

#Multiome human data
multiome_sce = readRDS(paste0(multiome_path, '/seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v5.rds'))
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))

multiome_sce


#Set up color scale
color_palette_1 = c('Astrocyte' = 'grey','Endo' = 'firebrick','Microglia' = 'darkred',
'Oligo' = 'darkgoldenrod','OPC' = 'cornsilk3')

color_palette_2 = c('Inhib.Thal' = 'dodgerblue','Excit.Thal' = 'indianred2','Thal' = 'lightsteelblue')

color_palette_3 = MetBrewer::met.brewer("Redon", n = 10)
names(color_palette_3) = c('LHb.1','LHb.1.3','LHb.1.3.4','LHb.2.7','LHb.4','LHb.7','MHb.1','MHb.1.2','MHb.2','MHb.3')

color_palette = c(color_palette_1, color_palette_2, color_palette_3)






#This doesn't have any saved dimension reductions
#Going back to when the data was converted from Seurat object to SCE, the reductions were not included
#See https://github.com/LieberInstitute/Hb_multiome/blob/master/code/08_spatial_registration_vs_multiome_snRNA-seq/01_multiome_rna_reference.R

#From whaT I can tell, the original umaps were generated in Seurat, using 30 nearest neighbores
#See https://github.com/LieberInstitute/Hb_multiome/blob/master/code/05_Clustering_ARCr/01_clustering_std_method.R

#Example, clust_knn was an argument passed to the batch job, assuming that matches the k30 in the data object name
#SeuratOBJ.1 <- RunUMAP(
#  SeuratOBJ.1,
#  n.neighbors = as.integer(clust_knn), # Default n.neighbors=30
#  nn.name = "weighted.nn",
#  reduction.name = "wnn.umap",
#  reduction.key = "wnnUMAP_"
#)


#Here, we're going to start with the default scran and scatter, runPCA, runUMAP functions
#The logcounts assay is the default used for PCA

#This is uncorrected data, across all the donors
# 1. Identify highly variable genes
dec <- modelGeneVar(multiome_sce)
hvg <- getTopHVGs(dec, n = 2000)

# 2. Run PCA on highly variable genes
multiome_sce <- runPCA(multiome_sce, subset_row = hvg)

# 3. Run UMAP on PCA space
multiome_sce <- runUMAP(multiome_sce, dimred = "PCA", n_dimred = 30)

# 4. Visualize
plotPCA(multiome_sce, colour_by = "mid_cluster")
plotUMAP(multiome_sce, colour_by = "mid_cluster")
plotUMAP(multiome_sce, colour_by = "orig.ident")

colnames(colData(multiome_sce))

table(multiome_sce$orig.ident, multiome_sce$mid_cluster)



#Exclude a donor, learn markers from the the remaining donors, project those markers onto the excluded donor.

#All the cell names are unique, so I can break apart by donor and then put it back together
table(duplicated(rownames(colData(multiome_sce))))

# Get unique donor IDs
donors <- unique(multiome_sce$orig.ident)

# Create a named list of SCE objects, one per donor
multiome_sce_list <- lapply(donors, function(donor) {
  multiome_sce[, multiome_sce$orig.ident == donor]
})
names(multiome_sce_list) <- donors

#Get the markers per donor
donor_markers_list = lapply(donors, function(donor) {
  MetaMarkers::compute_markers(assay(multiome_sce_list[[donor]], 'cpm'), multiome_sce_list[[donor]]$mid_cluster)
})
names(donor_markers_list) = donors


test_donor = "S08_Hb_r"

#Donor dimension reduction
test_sce = multiome_sce_list[[test_donor]]
dec <- modelGeneVar(test_sce)
hvg <- getTopHVGs(dec, n = 2000)

test_sce <- runPCA(test_sce, subset_row = hvg, exprs_values = "scaledata")
test_sce <- runUMAP(test_sce, dimred = "PCA", n_dimred = 20)



#Compute metamarkers excluding a donor
cross_donor_hab_markers = make_meta_markers(donor_markers_list[names(donor_markers_list) != test_donor], detailed_stats = TRUE)

#Project the metamarkers onto the excluded donor
top_markers = cross_donor_hab_markers %>% filter(cell_type != 'LHb.7' & rank <= 200)
ct_scores = score_cells(log1p(cpm(multiome_sce_list[[test_donor]])), top_markers)
ct_enrichment = compute_marker_enrichment(ct_scores)
ct_pred = assign_cells(ct_scores)

test_sce$cross_donor_pred_200 = ct_pred$predicted

top_markers = cross_donor_hab_markers %>% filter(cell_type != 'LHb.7' & rank <= 100)
ct_scores = score_cells(log1p(cpm(multiome_sce_list[[test_donor]])), top_markers)
ct_enrichment = compute_marker_enrichment(ct_scores)
ct_pred = assign_cells(ct_scores)

test_sce$cross_donor_pred_100 = ct_pred$predicted

top_markers = cross_donor_hab_markers %>% filter(cell_type != 'LHb.7' & rank <= 50) 
ct_scores = score_cells(log1p(cpm(multiome_sce_list[[test_donor]])), top_markers)
ct_enrichment = compute_marker_enrichment(ct_scores)
ct_pred = assign_cells(ct_scores)

test_sce$cross_donor_pred_50 = ct_pred$predicted


plotUMAP(test_sce, colour_by = "mid_cluster") +
  scale_color_manual(values = color_palette) + ggtitle('Original cluster annots')

#for (cluster in mid_res_annots) {
#  test_sce$highlight <- ifelse(test_sce$mid_cluster == cluster, cluster, "Other")
#  
#  p <- plotUMAP(test_sce, colour_by = "highlight") +
#    scale_color_manual(values = rlang::list2(!!cluster := color_palette[cluster], Other = "lightgrey")) +
#    ggtitle(cluster)
#  
#  print(p)
#}


plotUMAP(test_sce, colour_by = "cross_donor_pred_200") +
  scale_color_manual(values = color_palette) + ggtitle('Predicted cross-donor: top 200 markers')

plotUMAP(test_sce, colour_by = "cross_donor_pred_100") +
  scale_color_manual(values = color_palette) + ggtitle('Predicted cross-donor: top 100 markers')

plotUMAP(test_sce, colour_by = "cross_donor_pred_50") +
  scale_color_manual(values = color_palette) + ggtitle('Predicted cross-donor: top 50 markers')


table(test_sce$cross_donor_pred_100)



#Just going to go ahead with the cross-donor replicability clustering, just to compare to the integrated cluster labels

#Convert all the donor SCEs to seurats, do default seurat clustering, convert back to SCEs, and then metaneighbor across donors

#Convert first to seurats
multiome_seurat_list <- lapply(multiome_sce_list, function(sce) {
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
  multiome_seurat_list[[i]] <- RunPCA(multiome_seurat_list[[i]], features = VariableFeatures(object = multiome_seurat_list[[i]]))
  
  # Find neighbors and clusters
  multiome_seurat_list[[i]] <- FindNeighbors(multiome_seurat_list[[i]], dims = 1:20)
  multiome_seurat_list[[i]] <- FindClusters(multiome_seurat_list[[i]], resolution = .5)
  
  # UMAP
  multiome_seurat_list[[i]] <- RunUMAP(multiome_seurat_list[[i]], dims = 1:20)
  
  cat("Finished processing", names(multiome_seurat_list)[i], "\n")
}

#For trying different clustering resolutions
for (i in seq_along(multiome_seurat_list)) {
  multiome_seurat_list[[i]] <- FindClusters(multiome_seurat_list[[i]], resolution = .5)
}



#Check out the UMAPs for each donor
for (i in seq_along(multiome_seurat_list)) {
  p <- DimPlot(multiome_seurat_list[[i]], reduction = "umap", group.by = 'seurat_clusters', label = TRUE) + 
    ggtitle(names(multiome_seurat_list)[i])
  print(p)
  p <- DimPlot(multiome_seurat_list[[i]], reduction = "umap",group.by = 'mid_cluster' , label = TRUE) + 
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
mclusters = extractMetaClusters(MN_best_aurocs, threshold = .7)
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

logcounts(all_donor_sce) = log1p(cpm(all_donor_sce))
dec <- modelGeneVar(all_donor_sce)
hvg <- getTopHVGs(dec, n = 2000)

all_donor_sce <- runPCA(all_donor_sce, subset_row = hvg)
all_donor_sce <- runUMAP(all_donor_sce, dimred = "PCA", n_dimred = 20)

plotUMAP(all_donor_sce, colour_by = "mid_cluster") +
  scale_color_manual(values = color_palette) + ggtitle('Unintegrated mid_res clusters')

plotUMAP(all_donor_sce, colour_by = "meta_cluster") +
  ggtitle('Unintegrated meta clusters')








  