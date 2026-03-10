
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
color_palette_1 = c('Astrocyte' = 'grey','Endo' = 'firebrick','Microglia' = 'darkred',
'Oligo' = 'darkgoldenrod','OPC' = 'cornsilk3')

color_palette_2 = c('Inhib.Thal' = 'dodgerblue','Excit.Thal' = 'indianred2','Thal' = 'lightsteelblue')

color_palette_3 = MetBrewer::met.brewer("Redon", n = 10)
names(color_palette_3) = c('LHb.1','LHb.1.3','LHb.1.3.4','LHb.2.7','LHb.4','LHb.7','MHb.1','MHb.1.2','MHb.2','MHb.3')

color_palette = c(color_palette_1, color_palette_2, color_palette_3)


#Initial check on the co-expression of GABA and Glut markers across clusters
#Split by donor too

#Adapted from https://github.com/LieberInstitute/spatial_LS/blob/c68566780bc9a4e1aa3fe3f1c238d424bd4d8447/code/14_annotating_chromium_a-p/02_annotating.R#L150-L197

#Plotting co-expression of excitatory and inhibitory markers

pairwise_coexpression <- function(mat, genes, cluster_name) {
  # mat: genes x cells matrix for one cluster
  
  detected <- mat[genes, , drop = FALSE] > 0
  
  res <- expand.grid(gene1 = genes, gene2 = genes, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(percent = mean(detected[gene1, ] & detected[gene2, ]) * 100) %>%
    ungroup() %>%
    mutate(cluster = cluster_name)
  
  res
}

donors <- unique(colData(multiome_sce)[['orig.ident']])

for(donor in donors){
  sce_to_plot = multiome_sce[ , multiome_sce$orig.ident == donor]

  genes <- c('SLC32A1',"GAD1","GAD2","SLC17A6", "SLC17A7")
  expr_mat <- assay(sce_to_plot, "logcounts")

  cluster_to_annotate = "mid_cluster"
  clusters <- unique(colData(sce_to_plot)[[cluster_to_annotate]])

  coexp_df <- lapply(clusters, function(cl) {
    cells <- colData(sce_to_plot)[[cluster_to_annotate]] == cl
    mat_sub <- expr_mat[, cells, drop = FALSE]
    pairwise_coexpression(mat_sub, genes, cluster_name = cl)
  }) %>%
    bind_rows()

  coexp_df$gene1 <- factor(coexp_df$gene1, levels = genes)
  coexp_df$gene2 <- factor(coexp_df$gene2, levels = rev(genes))

  p <- ggplot(coexp_df, aes(x = gene1, y = gene2, fill = percent)) +
    geom_tile(color = "grey70", linewidth = 0.3) +
    facet_wrap(~ cluster, nrow = 5) +
    scale_fill_gradientn(
      colours = c("grey95","#f1e2c6", "#f1e2c6", "#f0c94a", "#df8b27", "#d92523", "#8b0d19"),
      values = c(0, 0.05, 0.20, 0.40, 0.60, 0.8, 1),
      limits = c(0, 100),
      breaks = c(0, 5, 20, 40, 60, 80, 100),
      name = "Percent of cells expressing\ntwo genes"
    ) +
    coord_equal() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    xlab(NULL) +
    ylab(NULL) + ggtitle(paste("Co-expression of GABA and Glut markers in", donor))

  print(p)
  
}


############################
#While the LHb.4 cluster does have some co-expression of GABA and Glut markers, they are super small percentages.
#I'm curious if these are just small subsets of thalamus cells that are present within the LHb.4 cluster
#Two pieces of evidence, the LHb.4 cluster is spatially the one closest to the thalamus
#And the GABA glut co-expression is the strongest in the thalamus clusters.

#So even a small amount of thalamus in the LHb.4 could account for the low levels of co GABA glut

##############################

#For the subset of GABA-Glut cells in a cluster, compare the expression of thalamus markers for those cells

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

#Make metamarkers, really just interested in the thalamus markers
cross_donor_hab_markers = make_meta_markers(donor_markers_list, detailed_stats = TRUE)



clust_filt = multiome_sce[ , multiome_sce$mid_cluster == 'LHb.4']
thal_filt = multiome_sce[ , multiome_sce$mid_cluster == 'Excit.Thal']
in_thal_filt = multiome_sce[ , multiome_sce$mid_cluster == 'Inhib.Thal']

#This gets the cells that co-express GAD2 and SLC17A6
gad1_expression = assay(clust_filt, 'cpm')['GAD1', ]
gad2_expression = assay(clust_filt, 'cpm')['GAD2', ]
vglut2_expression = assay(clust_filt, 'cpm')['SLC17A6', ]
vglut1_expression = assay(clust_filt, 'cpm')['SLC17A7', ]

coGabaGlut_index = gad2_expression * vglut2_expression > 0
mean(coGabaGlut_index)
sum(coGabaGlut_index)
#Get the average expression per cell of the top thalamus markers
thalamus_markers = cross_donor_hab_markers %>% filter(cell_type %in% c('Excit.Thal', 'Inhib.Thal') & rank <= 25) %>% pull(gene)
hab_markers = cross_donor_hab_markers %>% filter(cell_type %in% c('LHb.4') & rank <= 50) %>% pull(gene)

avg_thal_expr = colMeans(assay(clust_filt, 'cpm')[thalamus_markers, ])
avg_hab_expr = colMeans(assay(clust_filt, 'cpm')[hab_markers, ])
coGabaGlut_label = rep('Not', length = length(coGabaGlut_index))
coGabaGlut_label[coGabaGlut_index] = 'CoGABA-Glut'
tissue_label = rep('Habenula', length = length(coGabaGlut_index))
co_exp_df_1 = data.frame(avg_thal_expr = avg_thal_expr, avg_hab_expr = avg_hab_expr, 
  coGabaGlut_label = coGabaGlut_label, tissue = tissue_label)

avg_thal_expr = colMeans(assay(thal_filt, 'cpm')[thalamus_markers, ])
avg_hab_expr = colMeans(assay(thal_filt, 'cpm')[hab_markers, ])
coGabaGlut_label = rep('Excite.Thal', length = ncol(thal_filt))
tissue_label = rep('Thalamus', length = ncol(thal_filt))
co_exp_df_2 = data.frame(avg_thal_expr = avg_thal_expr, avg_hab_expr = avg_hab_expr,  
  coGabaGlut_label = coGabaGlut_label, tissue = tissue_label)

avg_thal_expr = colMeans(assay(in_thal_filt, 'cpm')[thalamus_markers, ])
avg_hab_expr = colMeans(assay(in_thal_filt, 'cpm')[hab_markers, ])
coGabaGlut_label = rep('Inhib.Thal', length = ncol(in_thal_filt))
tissue_label = rep('Thalamus', length = ncol(in_thal_filt))
co_exp_df_3 = data.frame(avg_thal_expr = avg_thal_expr, avg_hab_expr = avg_hab_expr,  
  coGabaGlut_label = coGabaGlut_label, tissue = tissue_label)

co_exp_df = rbind(co_exp_df_1, co_exp_df_2,co_exp_df_3)
co_exp_df$coGabaGlut_label = factor(co_exp_df$coGabaGlut_label, levels = c('Not', 'CoGABA-Glut', 'Excite.Thal', 'Inhib.Thal'))
ggplot(co_exp_df, aes(x = coGabaGlut_label, y = avg_thal_expr, fill = tissue)) +
  geom_violin(scale = 'width') +
  theme_bw() +
  xlab("Co-expression of GABA and Glut markers") +
  ylab("Average expression of top thalamus markers") +
  ggtitle("Thalamus marker expression in co-GABA-Glut cells vs others in LHb.4")

ggplot(co_exp_df, aes(x = coGabaGlut_label, y = avg_hab_expr, fill = tissue)) +
  geom_violin(scale = 'width') +
  theme_bw() +
  xlab("Co-expression of GABA and Glut markers") +
  ylab("Average expression of top LHb.4 markers") +
  ggtitle("Thalamus marker expression in co-GABA-Glut cells vs others in LHb.4")


############
#Okay so not obviously thalamus cells based on thalamus marker expression
###########








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








  