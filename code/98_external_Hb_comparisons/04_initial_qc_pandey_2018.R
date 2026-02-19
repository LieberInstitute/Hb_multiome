#Initial investigation into the public data for the zebrafish habenula data from the Pandey-Schier dataset
#Paper: https://www.sciencedirect.com/science/article/pii/S0960982218302252#abs0020


library(SingleCellExperiment)
library(Seurat)
library(MetaNeighbor)
library(dplyr)
library(ggplot2)
library(here)

here::here()

#Path to the Pandey 2018 data
pandey_10x_path = here('processed-data', '98_external_Hb_comparisons', 'Pandey_etal_2018_habenula_scseq', '10X')
list.files(pandey_10x_path)
pandey_ss_path = here('processed-data', '98_external_Hb_comparisons', 'Pandey_etal_2018_habenula_scseq', 'smartSeq')
list.files(pandey_ss_path)

#Path to save any generated data
new_data_path = here('processed-data', '98_external_Hb_comparisons', '04_initial_qc_pandey_2018')
#Path to plot directory
plot_path = here('plots', '98_external_Hb_comparisons', '04_initial_qc_pandey_2018')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


#The available data is in the from of count matrices
#Start with the 10X data, load up as seurat objects and do standard processing

# Read one of the count matrices
pandey_counts <- as.matrix(read.table(
  gzfile(paste0(pandey_10x_path, "/GSM2818521_larva_counts_matrix.txt.gz"))
))
# Create Seurat object
zeb_larva_hab_seu <- CreateSeuratObject(counts = pandey_counts, project = "Zebrafish_Habenula")
zeb_larva_hab_seu

#Adult 1 sample
pandey_counts <- as.matrix(read.table(
  gzfile(paste0(pandey_10x_path, "/GSM2818522_adultr1_counts_matrix.txt.gz"))
))
# Create Seurat object
zeb_adult1_hab_seu <- CreateSeuratObject(counts = pandey_counts, project = "Zebrafish_Habenula")
zeb_adult1_hab_seu


#Adult 2 sample
pandey_counts <- as.matrix(read.table(
  gzfile(paste0(pandey_10x_path, "/GSM2818523_adultr2_counts_matrix.txt.gz"))
))
# Create Seurat object
zeb_adult2_hab_seu <- CreateSeuratObject(counts = pandey_counts, project = "Zebrafish_Habenula")
zeb_adult2_hab_seu


#Not sure how mitochondrial genes are annotated in zebrafish
zeb_genes = rownames(zeb_larva_hab_seu)
mt_genes <- zeb_genes[grepl("^mt-|mt$|-mt-", zeb_genes, ignore.case = TRUE)]
#Looks like its still just MT-, so seurat default functions should work

zeb_larva_hab_seu[["percent.mt"]] <- PercentageFeatureSet(zeb_larva_hab_seu, pattern = "^MT-")
zeb_adult1_hab_seu[["percent.mt"]] <- PercentageFeatureSet(zeb_adult1_hab_seu, pattern = "^MT-")
zeb_adult2_hab_seu[["percent.mt"]] <- PercentageFeatureSet(zeb_adult2_hab_seu, pattern = "^MT-")

#Looks like these are actually pre-filtered, at least the adult samples match the 6% mt threshold reported in the paper
VlnPlot(zeb_larva_hab_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(zeb_adult1_hab_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(zeb_adult2_hab_seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Clustering steps for each individual donor

donor_list = list(zeb_larva_hab_seu, zeb_adult1_hab_seu, zeb_adult2_hab_seu)
names(donor_list) = c("Larva", "Adult1", "Adult2")
for (i in seq_along(donor_list)) {
  
  #CPM normalization
  donor_list[[i]] <- NormalizeData(donor_list[[i]], normalization.method = "RC", scale.factor = 1e6)

  # Find variable features
  donor_list[[i]] <- FindVariableFeatures(donor_list[[i]], selection.method = "vst", nfeatures = 2000)
  
  # Scale data
  all.genes <- rownames(donor_list[[i]])
  donor_list[[i]] <- ScaleData(donor_list[[i]], features = all.genes)
  
  # PCA
  donor_list[[i]] <- RunPCA(donor_list[[i]], features = VariableFeatures(object = donor_list[[i]]))
  
  # Find neighbors and clusters
  donor_list[[i]] <- FindNeighbors(donor_list[[i]], dims = 1:20)
  donor_list[[i]] <- FindClusters(donor_list[[i]], resolution = .8)
  
  # UMAP
  donor_list[[i]] <- RunUMAP(donor_list[[i]], dims = 1:20)
  
  cat("Finished processing", names(donor_list)[i], "\n")
}


#For trying different cluster resolutions
for (i in seq_along(donor_list)) {
  donor_list[[i]] <- FindClusters(donor_list[[i]], resolution = .8)
 
}


#Check out the UMAPs for each donor
#Generally seems to have the same structure, especially across the two adult samples
for (i in seq_along(donor_list)) {
  p <- DimPlot(donor_list[[i]], reduction = "umap", label = TRUE) + ggtitle(names(donor_list)[i])
  print(p)
}


#Save the donor-specific Seurat objects for future use
saveRDS(donor_list, file = paste0(new_data_path, '/individual_zebraf_seurat_objects_list.rds'))



#Convert to SingleCellExperiment objects for MetaNeighbor
donor_sce_list <- lapply(donor_list, as.SingleCellExperiment)
#Change the assays names to counts cpm scaledata
for (i in seq_along(donor_sce_list)) {
  names(assays(donor_sce_list[[i]])) <- c("counts", "cpm", "scaledata" )
}


#And run MetaNeighbor on the default Seurat clusters
#Get single SCE object
all_donor_sce = mergeSCE(donor_sce_list)
View(as.data.frame(colData(all_donor_sce)))

#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(dat = all_donor_sce, min_recurrence = 3, exp_labels = all_donor_sce$study_id)
length(global_hvgs)
keep_global_hvgs = global_hvgs[1:2000]


MN_aurocs = MetaNeighborUS(var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$seurat_clusters,
  fast_version = TRUE)

#Plot allby-all AUROC heatmap
plotHeatmap(MN_aurocs, 
  show_dendro = TRUE, 
  show_labels = TRUE, 
  cex = .5,
  title = "MetaNeighbor AUROCs for Pendey 2018 zebrafish")


#And the best versus next approach

MN_best_aurocs = MetaNeighborUS(var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$seurat_clusters,
  fast_version = TRUE,
one_vs_best = TRUE, symmetric_output = FALSE)

#Plot best_vs_next AUROC heatmap
plotHeatmap(MN_best_aurocs, 
  show_dendro = TRUE, 
  show_labels = TRUE, 
  cex = .5,
  title = "MetaNeighbor best_vs_next AUROCs for Pendey 2018 zebrafish")


#Get the metaclusters from the best vs next results, add those annotations to the full SCE object
mclusters = extractMetaClusters(MN_best_aurocs, threshold = .3)
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
# Check the result
head(all_donor_sce$meta_cluster)
table(names(all_donor_sce$meta_cluster), all_donor_sce$meta_cluster)


#Cluster the full SCE object and check out the larval vs adult annotations
#And the metacluster annotations

#Switch the full dataset to Seurat for the bubble plots
all_donor_seurat = as.Seurat(all_donor_sce, counts = "counts", data = "cpm")
all_donor_seurat

#Curious to see the meta-cluster annotations in the full dataset
#repeat the standard dim reduction from seurat

#HVGs
all_donor_seurat <- FindVariableFeatures(all_donor_seurat, selection.method = "vst", nfeatures = 2000)

#Scale data
all.genes <- rownames(all_donor_seurat)
all_donor_seurat <- ScaleData(all_donor_seurat, features = all.genes)

#PCA
all_donor_seurat  <- RunPCA(all_donor_seurat , features = VariableFeatures(object = all_donor_seurat ))
DimPlot(all_donor_seurat, reduction = "pca") + NoLegend()

#UMAP
all_donor_seurat  <- RunUMAP(all_donor_seurat , dims = 1:20)

p1 = DimPlot(all_donor_seurat , reduction = "umap", group.by = 'study_id', label = TRUE) + 
  ggtitle('Pendey 2018: batch annotations')
p1
ggsave(p1, filename = 'pendey_zebraf_hab_sample_batch_umap.pdf', path = plot_path,
device = 'pdf', width = 8, height = 7)

p2 = DimPlot(all_donor_seurat , reduction = "umap", group.by = 'meta_cluster', label = TRUE) + 
  ggtitle('Pendey 2018: MetaCluster annotations')
p2
ggsave(p2, filename = 'pendey_zebraf_hab_meta_cluster_umap.pdf', path = plot_path,
device = 'pdf', width = 8, height = 7)




#Save the combined SCE object with the donor-specific clusters as metadata for future use
saveRDS(all_donor_sce, file = paste0(new_data_path, '/all_donor_sce_with_denovo_clusters.rds'))
#all_donor_sce = readRDS(paste0(new_data_path, '/all_donor_sce_with_denovo_clusters.rds'))


#Save the MetaNeighbor results 
saveRDS(MN_aurocs, file = paste0(new_data_path, '/all_by_all_MN_aurocs.rds') )
saveRDS(MN_best_aurocs, file = paste0(new_data_path, '/best_vs_next_MN_aurocs.rds') )
#MN_best_aurocs = readRDS(paste0(new_data_path, '/best_vs_next_MN_aurocs.rds'))

paste0(plot_path, '/all_by_all_MN_aurocs_heatmap.pdf')
#Save the MetaNeighbor plots
pdf(paste0(plot_path, '/all_by_all_MN_aurocs_heatmap.pdf'), width = 10, height = 8)
#Plot allby-all AUROC heatmap
plotHeatmap(MN_aurocs, 
  show_dendro = TRUE, 
  show_labels = TRUE, 
  cex = .5,
  title = "MetaNeighbor AUROCs for Pendey 2018 zebrafish")

dev.off()


pdf(paste0(plot_path, '/best_vs_next_MN_aurocs_heatmap.pdf'), width = 10, height = 8)
#Plot best_vs_next AUROC heatmap
plotHeatmap(MN_best_aurocs, 
  show_dendro = TRUE, 
  show_labels = TRUE, 
  cex = .5,
  title = "MetaNeighbor best_vs_next AUROCs for Pendey 2018 zebrafish") 
dev.off()
















