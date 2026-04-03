#Now from the cross-species assessments, there is supporting evidence to possibly merge some of the multiome clusters
#As a recap on the cross-species clustering and potential cell-type labels

#Multiome clusters, mid resolution

#MHb.1 - Substance P
#MHb.2 - Cholinergic
#MHb.1.2 - unclear
#MHb.3 - small cluster that may be human specific, but will need to dig a bit deeper

#LHb.1 - Maps to the conserved ventral Lateral Hab cluster
#LHb.1.3 - Maps to the conserved ventral Lateral Hab cluster
#LHb.1.3.4 - Maps to the conserved ventral Lateral Hab cluster
#We would possibly merge these three

#LHb.4 - Maps to the conserved GABA/Glut Lateral Hab population
#LHb2.7 - Maps to the potentially mammal specific Lateral Hab cluster
#LHb.7 - Very small cluster mostly dominated by one donor, so probably need to drop


#From the ongoing discussions, some ideas on more evidence for merging or dropping clusters, 
#Pseudobulk PCA the clusters, see where they fall
#Maybe pseudobulk and PCA the human and mouse clusters together
#cell-type hierarchy based on transcriptomic similarity, see where the clusters fall
#A more extensive leave-one-out DE assessment and compare DE profiles. If excessivly similar, the maybe merge.


#Start with the pseudobulk PCA and cell-type hierarchy, simple and more straightforward.


library(SingleCellExperiment)
library(MetaMarkers)
library(qs2)
library(dplyr)
library(ggplot2)
library(scater)
library(scran)
library(here)

here::here()

#Path for new data generated
new_data_path = here('processed-data', '98_external_Hb_comparisons','11_multiome_cluster_merging')
#Path to plot directory
plot_path = here('plots', '98_external_Hb_comparisons', '11_multiome_cluster_merging')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '05_5_drop_doublets','01_drop_doublets_and_reDimReduce')

#Multiome human data
multiome_sce = qs_read(paste0(multiome_path, '/reprocessed_doubletRemoved_multiomeHab_SCE.qs2'))
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))

# Add in the putative inhibitory cluster annotations
inhib_data_path = here('processed-data', '99_donor_cluster_replicability','02_LHb4_investigation')
multiome_seurat_integrated = readRDS(paste0(inhib_data_path, '/multiome_LHb4_LHb7_integrated_seurat.rds'))
inhib_meta = multiome_seurat_integrated[[]]

inhib_barcodes_1 = rownames(inhib_meta)[inhib_meta$refined_mid_cluster == 'Putative_Inhib_LHb_4.1']
inhib_barcodes_2 = rownames(inhib_meta)[inhib_meta$refined_mid_cluster == 'Putative_Inhib_LHb_4.2']

# Add annotations
multiome_sce$refined_mid_cluster = multiome_sce$mid_cluster
multiome_sce$refined_mid_cluster[rownames(colData(multiome_sce)) %in% inhib_barcodes_1] = 'Putative_Inhib_LHb_4.1'
multiome_sce$refined_mid_cluster[rownames(colData(multiome_sce)) %in% inhib_barcodes_2] = 'Putative_Inhib_LHb_4.2'

table(multiome_sce$refined_mid_cluster, multiome_sce$mid_cluster)

#Check the donor distribtion, they have cells from all donors
table(multiome_sce$refined_mid_cluster, multiome_sce$orig.ident)

#Fine resolution annotations
multiome_sce$refined_cluster_ann = as.character(multiome_sce$cluster_ann)
multiome_sce$refined_cluster_ann[rownames(colData(multiome_sce)) %in% inhib_barcodes_1] = 'Putative_Inhib_LHb_4.1'
multiome_sce$refined_cluster_ann[rownames(colData(multiome_sce)) %in% inhib_barcodes_2] = 'Putative_Inhib_LHb_4.2'
table(multiome_sce$refined_cluster_ann )


#Pseudobulk the clusters from the counts, and then redo CPM on the summed counts
#Need the geneXcell (rowsXcolumns) matrix of counts
#Then need a matrix of cell annotations, cellsXannotations (rowsXcolumns)
#Matrix multiplication will get the sums of the counts for each gene for each cluster.

#Get the one-hot encoding of the cluster annotations for the multiome data, in this case the midcluster metadata column
cell_annot_matrix <- Matrix::sparse.model.matrix(~ 0 + refined_mid_cluster, data = colData(multiome_sce))
#Drops the 'mid_cluster' prefix
colnames(cell_annot_matrix) <- gsub("refined_mid_cluster", "", colnames(cell_annot_matrix))

#Get the pseudobulk counts
pseudobulk_counts = assay(multiome_sce, 'counts') %*% cell_annot_matrix

#Make a new SCE object with the pseudobulk counts
pseudobulk_coldata <- DataFrame(
  cell_type = colnames(pseudobulk_counts),
  n_cells = colSums(cell_annot_matrix)
)

pseudobulk_sce <- SingleCellExperiment(
  assays = list(counts = pseudobulk_counts),
  colData = pseudobulk_coldata
)

#Colors for plotting
color_palette = MetBrewer::met.brewer("Nizami", n = length(pseudobulk_sce$cell_type))
names(color_palette) = pseudobulk_sce$cell_type


#CPM the pseudobulk counts
assay(pseudobulk_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(pseudobulk_sce, 'counts'))

#Run PCA on the pseudobulked data
# 1. Log-transform the data 
assay(pseudobulk_sce, 'logcounts') <- log1p(assay(pseudobulk_sce, 'cpm'))

# 2. Identify highly variable genes
dec <- modelGeneVar(pseudobulk_sce)
hvg <- getTopHVGs(dec, n = 2000)

# 3. Run PCA
pseudobulk_sce <- runPCA(pseudobulk_sce, subset_row = hvg)

# 4. Visualize
plotPCA(pseudobulk_sce, colour_by = "cell_type", point_size = 5)

#For more control over plotting features
pca_data <- as.data.frame(reducedDim(pseudobulk_sce, "PCA")[, 1:2])
pca_data$cell_type <- pseudobulk_sce$cell_type

# Get percent variance explained
pca_attr <- attr(reducedDim(pseudobulk_sce, "PCA"), "percentVar")
pc1_var <- round(pca_attr[1], 1)
pc2_var <- round(pca_attr[2], 1)

p_pca_all_multiome = ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_type, label = cell_type)) +
  geom_point(size = 5) +
  ggrepel::geom_label_repel(show.legend = FALSE,
                            box.padding = 0.5, 
                            point.padding = 0.5,
                            max.overlaps = Inf,
                            min.segment.length = 0) +
  scale_color_manual(values = color_palette, name = 'Cell-type') +
  theme_bw() +
  labs(x = paste0("PC1 (", pc1_var, "%)"),
       y = paste0("PC2 (", pc2_var, "%)"),
       title = "Pseudobulk PCA of all Multiome Clusters")

p_pca_all_multiome

#And now just look at the neurons
pseudobulk_neurons_sce = pseudobulk_sce[, grepl('^MHb|^LHb|^Putative', pseudobulk_sce$cell_type)]

dec_neurons <- modelGeneVar(pseudobulk_neurons_sce)
hvg_neurons <- getTopHVGs(dec_neurons, n = 2000)

pseudobulk_neurons_sce <- runPCA(pseudobulk_neurons_sce, subset_row = hvg_neurons)

#Visualize
pca_data_neurons <- as.data.frame(reducedDim(pseudobulk_neurons_sce, "PCA")[, 1:2])
pca_data_neurons$cell_type <- pseudobulk_neurons_sce$cell_type

# Get percent variance explained
pca_attr_neurons <- attr(reducedDim(pseudobulk_neurons_sce, "PCA"), "percentVar")
pc1_var_neurons <- round(pca_attr_neurons[1], 1)
pc2_var_neurons <- round(pca_attr_neurons[2], 1)

p_pca_neuron_multiome = ggplot(pca_data_neurons, aes(x = PC1, y = PC2, color = cell_type, label = cell_type)) +
  geom_point(size = 5) +
  ggrepel::geom_label_repel(show.legend = FALSE,
                            box.padding = 0.5, 
                            point.padding = 0.5,
                            max.overlaps = Inf,
                            min.segment.length = 0) +
  scale_color_manual(values = color_palette, name = 'Cell-type') +
  theme_bw() +
  labs(x = paste0("PC1 (", pc1_var_neurons, "%)"),
       y = paste0("PC2 (", pc2_var_neurons, "%)"),
       title = "Pseudobulk PCA of neuronal Multiome Clusters")

p_pca_neuron_multiome 

ggsave(p_pca_all_multiome, filename = 'multiome_pseudobulk_all_cluster_pca.pdf', path = plot_path, useDingbats = F,
width = 8, height = 6, device = 'pdf')

ggsave(p_pca_neuron_multiome, filename = 'multiome_pseudobulk_neuron_cluster_pca.pdf', path = plot_path, useDingbats = F,
width = 8, height = 6, device = 'pdf')



#And what about a pseudobulk pca with the fine resolution clusters
#Get the one-hot encoding of the cluster annotations for the multiome data, in this case the midcluster metadata column
cell_annot_matrix <- Matrix::sparse.model.matrix(~ 0 + refined_cluster_ann, data = colData(multiome_sce))
colnames(cell_annot_matrix) <- gsub("refined_cluster_ann", "", colnames(cell_annot_matrix))

#Get the pseudobulk counts
pseudobulk_counts = assay(multiome_sce, 'counts') %*% cell_annot_matrix

#Make a new SCE object with the pseudobulk counts
pseudobulk_coldata <- DataFrame(
  cell_type = colnames(pseudobulk_counts),
  n_cells = colSums(cell_annot_matrix)
)

pseudobulk_sce <- SingleCellExperiment(
  assays = list(counts = pseudobulk_counts),
  colData = pseudobulk_coldata
)

#Colors for plotting
color_palette = MetBrewer::met.brewer("Nizami", n = length(pseudobulk_sce$cell_type))
names(color_palette) = pseudobulk_sce$cell_type


#CPM the pseudobulk counts
assay(pseudobulk_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(pseudobulk_sce, 'counts'))

#Run PCA on the pseudobulked data
# 1. Log-transform the data 
assay(pseudobulk_sce, 'logcounts') <- log1p(assay(pseudobulk_sce, 'cpm'))

# 2. Identify highly variable genes
dec <- modelGeneVar(pseudobulk_sce)
hvg <- getTopHVGs(dec, n = 2000)

# 3. Run PCA
pseudobulk_sce <- runPCA(pseudobulk_sce, subset_row = hvg)

# 4. Visualize
plotPCA(pseudobulk_sce, colour_by = "cell_type", point_size = 5)

#For more control over plotting features
pca_data <- as.data.frame(reducedDim(pseudobulk_sce, "PCA")[, 1:2])
pca_data$cell_type <- pseudobulk_sce$cell_type

# Get percent variance explained
pca_attr <- attr(reducedDim(pseudobulk_sce, "PCA"), "percentVar")
pc1_var <- round(pca_attr[1], 1)
pc2_var <- round(pca_attr[2], 1)

p_pca_all_multiome_fine = ggplot(pca_data, aes(x = PC1, y = PC2, color = cell_type, label = cell_type)) +
  geom_point(size = 5) +
  ggrepel::geom_label_repel(show.legend = FALSE,
                            box.padding = 0.5, 
                            point.padding = 0.5,
                            max.overlaps = Inf,
                            min.segment.length = 0) +
  scale_color_manual(values = color_palette, name = 'Cell-type') +
  theme_bw() +
  labs(x = paste0("PC1 (", pc1_var, "%)"),
       y = paste0("PC2 (", pc2_var, "%)"),
       title = "Pseudobulk PCA of all Multiome Fine Clusters")

p_pca_all_multiome_fine

ggsave(p_pca_all_multiome_fine, filename = 'multiome_pseudobulk_all_fine_cluster_pca.pdf', path = plot_path, useDingbats = F,
width = 12, height = 10, device = 'pdf')




#And now build the cluster hierarchy 
#Compute markers for each annotation, using the full dataset
#take the top 200 markers for each cluster, get the unique subset of these genes
#Compute expression centroids, correlations across those centroids, and then hierarchical clustering on those correlations
top_markers = compute_markers(assay(multiome_sce, 'cpm'), multiome_sce$refined_mid_cluster)

top_markers %>% group_by(cell_type) %>% slice_max(auroc, n= 10) %>% View()


#Unique list of the top markers
top_marker_sub = top_markers %>% group_by(cell_type) %>% slice_max(auroc, n= 500) %>%
  pull(gene) %>% unique()

#Expression data for those genes
gene_filt = rownames(multiome_sce) %in% top_marker_sub
exp_data = assay(multiome_sce, 'cpm')[gene_filt, ]
#Compute centroids of gene expression per subclass
exp_data = as.data.frame(t(as.matrix(exp_data)))
exp_data = exp_data %>% mutate(celltype = multiome_sce$refined_mid_cluster)
centroids = exp_data %>% group_by(celltype) %>% summarize(across(which(colnames(exp_data)!= 'celltype'), median)) %>% as.data.frame()
#Get back into genes on rows and subclass on columns
rownames(centroids) = centroids$celltype
centroids = t(centroids[ ,2:ncol(centroids)])

#Compute a distance matrix from the correlations of the centroids
centroid_corr = cor(centroids, method = 'spearman')

#Hierarchical clusting on the correlation distance matrix
hclust_avg <- hclust(as.dist(1-centroid_corr), method = 'average')

#Visualize the correlation matrix the dendrogram is derived from

viridis_map = circlize::colorRamp2(seq(0, 1, length.out = 100),
                                  viridis::rocket(100))

cent_hm_500 = ComplexHeatmap::Heatmap(centroid_corr, col = viridis_map, name = 'spearman' , show_row_dend = FALSE,
                                  clustering_distance_columns = function(m) dist(1-m), clustering_method_columns = "average",
                                  clustering_distance_rows = function(m) dist(1-m), clustering_method_rows = "average",
                                  column_dend_height = unit(3, "cm"),
                                  column_title = "Multiome cluster taxonomy: top 500 markers")
cent_hm_500 = ComplexHeatmap::draw(cent_hm_500)


#Pull out all the pairwise points and highlight where the LHb1 clusters are in the distribution
# Get upper triangle indices (excluding diagonal)
upper_tri <- which(upper.tri(centroid_corr), arr.ind = TRUE)

# Create the pairwise comparison dataframe
pairwise_corr <- data.frame(
  row_name = rownames(centroid_corr)[upper_tri[, 1]],
  col_name = colnames(centroid_corr)[upper_tri[, 2]],
  correlation = centroid_corr[upper_tri]
) %>%
  mutate(comparison = paste(row_name, col_name, sep = " vs "))

# Flag comparisons involving only LHb.1, LHb.1.3, and LHb.1.3.4
lhb_clusters <- c("LHb.1", "LHb.1.3", "LHb.1.3.4")
pairwise_corr <- pairwise_corr %>%
  mutate(label = ifelse(row_name %in% lhb_clusters & col_name %in% lhb_clusters, comparison, NA))

pairwise_corr_500 = pairwise_corr

# Add jittered x position to the data
set.seed(123)  # For reproducibility
pairwise_corr$x_jitter <- jitter(rep(1, nrow(pairwise_corr)), amount = 0.05)

# Create boxplot with points
ggplot(pairwise_corr, aes(x = x_jitter, y = correlation)) +
  geom_boxplot(aes(x = 1), width = .3) +
  geom_point(alpha = 0.5, size = 2) +
  ggrepel::geom_label_repel(aes(label = label), 
                   na.rm = TRUE,
                   box.padding = 0.5,
                   point.padding = 0.5,
                   min.segment.length = 0) +
  theme_bw() +
  scale_x_continuous(breaks = NULL) +
  labs(x = "", y = "Spearman Correlation", 
       title = "Pairwise cluster correlations")


#And the top 200 markers

#Unique list of the top markers
top_marker_sub = top_markers %>% group_by(cell_type) %>% slice_max(auroc, n= 200) %>%
  pull(gene) %>% unique()

#Expression data for those genes
gene_filt = rownames(multiome_sce) %in% top_marker_sub
exp_data = assay(multiome_sce, 'cpm')[gene_filt, ]
#Compute centroids of gene expression per subclass
exp_data = as.data.frame(t(as.matrix(exp_data)))
exp_data = exp_data %>% mutate(celltype = multiome_sce$refined_mid_cluster)
centroids = exp_data %>% group_by(celltype) %>% summarize(across(which(colnames(exp_data)!= 'celltype'), median)) %>% as.data.frame()
#Get back into genes on rows and subclass on columns
rownames(centroids) = centroids$celltype
centroids = t(centroids[ ,2:ncol(centroids)])

#Compute a distance matrix from the correlations of the centroids
centroid_corr = cor(centroids, method = 'spearman')

#Hierarchical clusting on the correlation distance matrix
hclust_avg <- hclust(as.dist(1-centroid_corr), method = 'average')

#Visualize the correlation matrix the dendrogram is derived from

viridis_map = circlize::colorRamp2(seq(0, 1, length.out = 100),
                                  viridis::rocket(100))

cent_hm_200 = ComplexHeatmap::Heatmap(centroid_corr, col = viridis_map, name = 'spearman' , show_row_dend = FALSE,
                                  clustering_distance_columns = function(m) dist(1-m), clustering_method_columns = "average",
                                  clustering_distance_rows = function(m) dist(1-m), clustering_method_rows = "average",
                                  column_dend_height = unit(3, "cm"),
                                  column_title = "Multiome cluster taxonomy: top 200 markers")
cent_hm_200 = ComplexHeatmap::draw(cent_hm_200)


#And the top 50 markers
#Unique list of the top markers
top_marker_sub = top_markers %>% group_by(cell_type) %>% slice_max(auroc, n= 50) %>%
  pull(gene) %>% unique()

#Expression data for those genes
gene_filt = rownames(multiome_sce) %in% top_marker_sub
exp_data = assay(multiome_sce, 'cpm')[gene_filt, ]
#Compute centroids of gene expression per subclass
exp_data = as.data.frame(t(as.matrix(exp_data)))
exp_data = exp_data %>% mutate(celltype = multiome_sce$refined_mid_cluster)
centroids = exp_data %>% group_by(celltype) %>% summarize(across(which(colnames(exp_data)!= 'celltype'), median)) %>% as.data.frame()
#Get back into genes on rows and subclass on columns
rownames(centroids) = centroids$celltype
centroids = t(centroids[ ,2:ncol(centroids)])

#Compute a distance matrix from the correlations of the centroids
centroid_corr = cor(centroids, method = 'spearman')

#Hierarchical clusting on the correlation distance matrix
hclust_avg <- hclust(as.dist(1-centroid_corr), method = 'average')

#Visualize the correlation matrix the dendrogram is derived from

viridis_map = circlize::colorRamp2(seq(0, 1, length.out = 100),
                                  viridis::rocket(100))

cent_hm_50 = ComplexHeatmap::Heatmap(centroid_corr, col = viridis_map, name = 'spearman' , show_row_dend = FALSE,
                                  clustering_distance_columns = function(m) dist(1-m), clustering_method_columns = "average",
                                  clustering_distance_rows = function(m) dist(1-m), clustering_method_rows = "average",
                                  column_dend_height = unit(3, "cm"),
                                  column_title = "Multiome cluster taxonomy: top 50 markers")
cent_hm_50 = ComplexHeatmap::draw(cent_hm_50)


#And the top 25 markers
#Unique list of the top markers
top_marker_sub = top_markers %>% group_by(cell_type) %>% slice_max(auroc, n= 25) %>%
  pull(gene) %>% unique()

#Expression data for those genes
gene_filt = rownames(multiome_sce) %in% top_marker_sub
exp_data = assay(multiome_sce, 'cpm')[gene_filt, ]
#Compute centroids of gene expression per subclass
exp_data = as.data.frame(t(as.matrix(exp_data)))
exp_data = exp_data %>% mutate(celltype = multiome_sce$refined_mid_cluster)
centroids = exp_data %>% group_by(celltype) %>% summarize(across(which(colnames(exp_data)!= 'celltype'), median)) %>% as.data.frame()
#Get back into genes on rows and subclass on columns
rownames(centroids) = centroids$celltype
centroids = t(centroids[ ,2:ncol(centroids)])

#Compute a distance matrix from the correlations of the centroids
centroid_corr = cor(centroids, method = 'spearman')

#Hierarchical clusting on the correlation distance matrix
hclust_avg <- hclust(as.dist(1-centroid_corr), method = 'average')

#Visualize the correlation matrix the dendrogram is derived from

viridis_map = circlize::colorRamp2(seq(0, 1, length.out = 100),
                                  viridis::rocket(100))

cent_hm_25 = ComplexHeatmap::Heatmap(centroid_corr, col = viridis_map, name = 'spearman' , show_row_dend = FALSE,
                                  clustering_distance_columns = function(m) dist(1-m), clustering_method_columns = "average",
                                  clustering_distance_rows = function(m) dist(1-m), clustering_method_rows = "average",
                                  column_dend_height = unit(3, "cm"),
                                  column_title = "Multiome cluster taxonomy: top 25 markers")
cent_hm_25 = ComplexHeatmap::draw(cent_hm_25)


#Pull out all the pairwise points and highlight where the LHb1 clusters are in the distribution
# Get upper triangle indices (excluding diagonal)
upper_tri <- which(upper.tri(centroid_corr), arr.ind = TRUE)

# Create the pairwise comparison dataframe
pairwise_corr <- data.frame(
  row_name = rownames(centroid_corr)[upper_tri[, 1]],
  col_name = colnames(centroid_corr)[upper_tri[, 2]],
  correlation = centroid_corr[upper_tri]
) %>%
  mutate(comparison = paste(row_name, col_name, sep = " vs "))

# Flag comparisons involving only LHb.1, LHb.1.3, and LHb.1.3.4
lhb_clusters <- c("LHb.1", "LHb.1.3", "LHb.1.3.4")
pairwise_corr <- pairwise_corr %>%
  mutate(label = ifelse(row_name %in% lhb_clusters & col_name %in% lhb_clusters, comparison, NA))

pairwise_corr_200 = pairwise_corr

# Add jittered x position to the data
set.seed(123)  # For reproducibility
pairwise_corr$x_jitter <- jitter(rep(1, nrow(pairwise_corr)), amount = 0.05)

# Create boxplot with points
ggplot(pairwise_corr, aes(x = x_jitter, y = correlation)) +
  geom_boxplot(aes(x = 1), width = .3) +
  geom_point(alpha = 0.5, size = 2) +
  ggrepel::geom_label_repel(aes(label = label), 
                   na.rm = TRUE,
                   box.padding = 0.5,
                   point.padding = 0.5,
                   min.segment.length = 0) +
  theme_bw() +
  scale_x_continuous(breaks = NULL) +
  labs(x = "", y = "Spearman Correlation", 
       title = "Pairwise cluster correlations")


pairwise_corr_200 %>% arrange(desc(correlation)) %>% head(10)
pairwise_corr_500 %>% arrange(desc(correlation)) %>% head(10)


pdf(file = paste0(plot_path, '/multiome_cluster_taxonomy_top25_markers.pdf'), useDingbats = F, 
width = 6, height = 6)
cent_hm_25 
dev.off()

pdf(file = paste0(plot_path, '/multiome_cluster_taxonomy_top50_markers.pdf'), useDingbats = F, 
width = 6, height = 6)
cent_hm_50 
dev.off()

pdf(file = paste0(plot_path, '/multiome_cluster_taxonomy_top200_markers.pdf'), useDingbats = F, 
width = 6, height = 6)
cent_hm_200 
dev.off()

pdf(file = paste0(plot_path, '/multiome_cluster_taxonomy_top500_markers.pdf'), useDingbats = F, 
width = 6, height = 6)
cent_hm_500 
dev.off()



#And what about a fine resolution hierarcy, just try 50 markers

top_markers = compute_markers(assay(multiome_sce, 'cpm'), multiome_sce$refined_cluster_ann)

#Unique list of the top markers
top_marker_sub = top_markers %>% group_by(cell_type) %>% slice_max(auroc, n= 50) %>%
  pull(gene) %>% unique()

#Expression data for those genes
gene_filt = rownames(multiome_sce) %in% top_marker_sub
exp_data = assay(multiome_sce, 'cpm')[gene_filt, ]
#Compute centroids of gene expression per subclass
exp_data = as.data.frame(t(as.matrix(exp_data)))
exp_data = exp_data %>% mutate(celltype = multiome_sce$refined_cluster_ann)
centroids = exp_data %>% group_by(celltype) %>% summarize(across(which(colnames(exp_data)!= 'celltype'), median)) %>% as.data.frame()
#Get back into genes on rows and subclass on columns
rownames(centroids) = centroids$celltype
centroids = t(centroids[ ,2:ncol(centroids)])

#Compute a distance matrix from the correlations of the centroids
centroid_corr = cor(centroids, method = 'spearman')

#Hierarchical clusting on the correlation distance matrix
hclust_avg <- hclust(as.dist(1-centroid_corr), method = 'average')

#Visualize the correlation matrix the dendrogram is derived from

viridis_map = circlize::colorRamp2(seq(0, 1, length.out = 100),
                                  viridis::rocket(100))

cent_fine_hm_50 = ComplexHeatmap::Heatmap(centroid_corr, col = viridis_map, name = 'spearman' , show_row_dend = FALSE,
                                  clustering_distance_columns = function(m) dist(1-m), clustering_method_columns = "average",
                                  clustering_distance_rows = function(m) dist(1-m), clustering_method_rows = "average",
                                  column_dend_height = unit(3, "cm"),
                                  column_title = "Multiome fine cluster taxonomy: top 50 markers")
cent_fine_hm_50 = ComplexHeatmap::draw(cent_fine_hm_50)

pdf(file = paste0(plot_path, '/multiome_fine_cluster_taxonomy_top50_markers.pdf'), useDingbats = F, 
width = 6, height = 6)
cent_fine_hm_50 
dev.off()





