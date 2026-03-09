
#One of the main conclusions from the cross-species comparisons, is that the current MHb.3 population may be human specific
#It is a small cluster, but almost equally contributed to across donors, and does not have a strong match to mouse or zebrafish clusters.

#This analysis step will do a bit more digging into the MHb.3 cluster.
#Central idea, is to see if this cluster arises obviously when clustering at a per-donor level.


library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
library(dplyr)
library(ggplot2)
library(scater)
library(scran)
library(here)

here::here()

#Path for new data generated
new_data_path = here('processed-data', '99_donor_cluster_replicability','01_per_donor_clustering')
#Path to plot directory
plot_path = here('plots', '99_donor_cluster_replicability','01_per_donor_clustering')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '08_spatial_registration_vs_multiome_snRNA-seq','mid')

#Multiome human data
multiome_sce = readRDS(paste0(multiome_path, '/seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v5.rds'))
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))

multiome_sce

#Set up color scale
mid_res_annots = names(table(multiome_sce$mid_cluster))
#color_palette_1 = MetBrewer::met.brewer("Monet", n = 5)
#names(color_palette_1) = c('Astrocyte','Endo','Microglia','Oligo','OPC')

color_palette_1 = c('Astrocyte' = 'grey','Endo' = 'firebrick','Microglia' = 'darkred',
'Oligo' = 'darkgoldenrod','OPC' = 'cornsilk3')

#color_palette_2 = MetBrewer::met.brewer("Ingres", n = 3)
#names(color_palette_2) = c('Inhib.Thal','Excit.Thal','Thal')

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



#And what does seurat dimreduc look like?

test_seurat = as.Seurat(test_sce, counts = 'counts', data = 'cpm')
test_seurat  <- FindVariableFeatures(test_seurat , selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(test_seurat )
test_seurat  <- ScaleData(test_seurat , features = all.genes)

test_seurat <- RunPCA(test_seurat, features = VariableFeatures(object = test_seurat))

test_seurat <- FindNeighbors(test_seurat, dims = 1:20)
test_seurat <- FindClusters(test_seurat, resolution = .5)

test_seurat <- RunUMAP(test_seurat, dims = 1:20)

DimPlot(test_seurat, reduction = "umap", group.by = "mid_cluster", cols = color_palette, pt.size = 2) + 
  ggtitle('Seurat UMAP: original cluster annots')
DimPlot(test_seurat, reduction = "umap", group.by = 'seurat_clusters', pt.size = 2) + 
  ggtitle('Seurat UMAP')

test_seurat@meta.data %>% View()






