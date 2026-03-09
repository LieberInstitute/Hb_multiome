
#One of the main conclusions from the cross-species comparisons, is that the current MHb.3 population may be human specific
#It is a small cluster, but almost equally contributed to across donors, and does not have a strong match to mouse or zebrafish clusters.

#This analysis step will do a bit more digging into the MHb.3 cluster.
#Central idea, is to see if this cluster arises obviously when clustering at a per-donor level.


library(SingleCellExperiment)
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
color_palette = MetBrewer::met.brewer("Redon", n = length(mid_res_annots))
names(color_palette) = mid_res_annots


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

table(multiome_sce$orig.ident)



#What about a single donor?
#Donor S11_Hb_r has the most cells

s11_sce = multiome_sce[, multiome_sce$orig.ident == 'S11_Hb_r']

dec <- modelGeneVar(s11_sce)
hvg <- getTopHVGs(dec, n = 2000)

s11_sce <- runPCA(s11_sce, subset_row = hvg)

s11_sce <- runUMAP(s11_sce, dimred = "PCA", n_dimred = 30)

# 4. Visualize
plotPCA(s11_sce, colour_by = "mid_cluster")
plotUMAP(s11_sce, colour_by = "mid_cluster") +
  scale_color_manual(values = color_palette)

for (cluster in mid_res_annots) {
  s11_sce$highlight <- ifelse(s11_sce$mid_cluster == cluster, cluster, "Other")
  
  p <- plotUMAP(s11_sce, colour_by = "highlight") +
    scale_color_manual(values = rlang::list2(!!cluster := color_palette[cluster], Other = "lightgrey")) +
    ggtitle(cluster)
  
  print(p)
}
















