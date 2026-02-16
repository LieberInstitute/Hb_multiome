#Script for initial basic QC and attempting default Seurat clustering of the Wallace 2019 mouse data.

library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(MetaNeighbor)

#Load the data
#Cell count matches the reported 7506 final cells used in the paper
hab_batch1_data = readRDS('processed-data/98_external_Hb_comparisons/Wallace_etal_2019_habenula_scseq/hab_batch1.rds')
hab_batch1_data

#Made with Seurat v2.3.4
hab_batch1_data@version

#Update the object to v3
hab_batch1_data <- UpdateSeuratObject(hab_batch1_data)


#Sample metadata is likely hidden in the cell barcodes, though never explained by the authors
barcodes = rownames(hab_batch1_data@meta.data)
barcode_df <- as.data.frame(do.call(rbind, strsplit(barcodes, split = '_', fixed = TRUE)))
head(barcode_df)

#Looks like V2 is donor and V4 is left or right habenula
#Seems like there are only 4 donors, but the paper reports 6. The cell number matches though, still 7506
barcode_df %>% group_by(V1, V2, V4) %>% 
  summarise(count = n())
#Add the donor and left/right info to the metadata

hab_batch1_data$putative_donor = barcode_df$V2
hab_batch1_data$putative_LR = barcode_df$V4





#Follow along basic Seurat cluster workflow, as per https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
VlnPlot(hab_batch1_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

#Looks like already pre-QC'd data.
#Mitochondrial % between 0-10
#counts between 501 and 17787
#features between 201 and 5276
min(hab_batch1_data$percent.mito)
min(hab_batch1_data$nCount_RNA)
min(hab_batch1_data$nFeature_RNA)

max(hab_batch1_data$percent.mito)
max(hab_batch1_data$nCount_RNA)
max(hab_batch1_data$nFeature_RNA)

#Normalize the data, sticking with CPM for now 
hab_batch1_data <- NormalizeData(hab_batch1_data, normalization.method = "RC", scale.factor = 1e6)
hab_batch1_data@assays$RNA@data[1:5,1:5]

#HVGs
hab_batch1_data <- FindVariableFeatures(hab_batch1_data, selection.method = "vst", nfeatures = 2000)

#Scale data
all.genes <- rownames(hab_batch1_data)
hab_batch1_data <- ScaleData(hab_batch1_data, features = all.genes)

#PCA
hab_batch1_data  <- RunPCA(hab_batch1_data , features = VariableFeatures(object = hab_batch1_data ))
DimPlot(hab_batch1_data, reduction = "pca") + NoLegend()
ElbowPlot(hab_batch1_data, ndims = 50)


#Clustering
hab_batch1_data  <- FindNeighbors(hab_batch1_data , dims = 1:20)
#Start with default resolution, which is 0.8
hab_batch1_data  <- FindClusters(hab_batch1_data )


#UMAP
hab_batch1_data  <- RunUMAP(hab_batch1_data , dims = 1:20)
DimPlot(hab_batch1_data , reduction = "umap")

#Check out the sample metadata to see if it matches to clear batch effects
DimPlot(hab_batch1_data, group.by = "putative_donor")
DimPlot(hab_batch1_data, group.by = "putative_LR")

#Yup, as I thought. Obvious batch structure between donors 161103 + 161105 and donors 160822 + 161102
#The L and R habenula are well mixed.

#One simple analysis to try, just as practice to get things up and running, is a per-donor cluster plus MetaNeighbor assessment



#Break up by donor and then cluster each separately
donors <- unique(hab_batch1_data$putative_donor)
donor_list <- lapply(donors, function(donor) {
  subset(hab_batch1_data, subset = putative_donor == donor)
})
names(donor_list) <- paste0("hab_", donors)



#Clustering steps for each individual donor

for (i in seq_along(donor_list)) {
  # Find variable features
  donor_list[[i]] <- FindVariableFeatures(donor_list[[i]], selection.method = "vst", nfeatures = 2000)
  
  # Scale data
  all.genes <- rownames(donor_list[[i]])
  donor_list[[i]] <- ScaleData(donor_list[[i]], features = all.genes)
  
  # PCA
  donor_list[[i]] <- RunPCA(donor_list[[i]], features = VariableFeatures(object = donor_list[[i]]))
  
  # Find neighbors and clusters
  donor_list[[i]] <- FindNeighbors(donor_list[[i]], dims = 1:20)
  donor_list[[i]] <- FindClusters(donor_list[[i]])
  
  # UMAP
  donor_list[[i]] <- RunUMAP(donor_list[[i]], dims = 1:20)
  
  cat("Finished processing", names(donor_list)[i], "\n")
}

#Check out the UMAPs for each donor
for (i in seq_along(donor_list)) {
  p <- DimPlot(donor_list[[i]], reduction = "umap", label = TRUE) + ggtitle(names(donor_list)[i])
  print(p)
}


#Save the donor-specific Seurat objects for future use
saveRDS(donor_list, file = 'processed-data/98_external_Hb_comparisons/02_qc_and_clust_wallace_2019/individual_donor_seurat_objects_list.rds')




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
global_hvgs = variableGenes(dat = all_donor_sce, min_recurrence = 4, exp_labels = all_donor_sce$study_id)
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
  title = "MetaNeighbor AUROCs for Wallace 2019 donors")

#Get the reciprocal best hits
wallace_MN_top_hits_df = topHits(MN_aurocs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$seurat_clusters,
  threshold = .9)
View(wallace_MN_top_hits_df)


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
  title = "MetaNeighbor best_vs_next AUROCs for Wallace 2019 donors")




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
# Check the result
head(all_donor_sce$meta_cluster)
table(names(all_donor_sce$meta_cluster), all_donor_sce$meta_cluster)


#Save the combined SCE object with the donor-specific clusters as metadata for future use
saveRDS(all_donor_sce, file = 'processed-data/98_external_Hb_comparisons/02_qc_and_clust_wallace_2019/all_donor_sce_with_denovo_clusters.rds')
#all_donor_sce = readRDS('processed-data/98_external_Hb_comparisons/02_qc_and_clust_wallace_2019/all_donor_sce_with_denovo_clusters.rds')


#Save the MetaNeighbor results 
saveRDS(MN_aurocs, file = 'processed-data/98_external_Hb_comparisons/02_qc_and_clust_wallace_2019/all_by_all_MN_aurocs.rds' )
saveRDS(MN_best_aurocs, file = 'processed-data/98_external_Hb_comparisons/02_qc_and_clust_wallace_2019/best_vs_next_MN_aurocs.rds' )
#MN_best_aurocs = readRDS('processed-data/98_external_Hb_comparisons/02_qc_and_clust_wallace_2019/best_vs_next_MN_aurocs.rds' )


#Save the MetaNeighbor plots
pdf('plots/98_external_Hb_comparisons/02_qc_and_clust_wallace_2019/all_by_all_MN_aurocs_heatmap.pdf', width = 10, height = 8)
#Plot allby-all AUROC heatmap
plotHeatmap(MN_aurocs, 
  show_dendro = TRUE, 
  show_labels = TRUE, 
  cex = .5,
  title = "MetaNeighbor AUROCs for Wallace 2019 donors")

dev.off()

pdf('plots/98_external_Hb_comparisons/02_qc_and_clust_wallace_2019/best_vs_next_MN_aurocs_heatmap.pdf', width = 10, height = 8)
#Plot best_vs_next AUROC heatmap
plotHeatmap(MN_best_aurocs, 
  show_dendro = TRUE, 
  show_labels = TRUE, 
  cex = .5,
  title = "MetaNeighbor best_vs_next AUROCs for Wallace 2019 donors") 
dev.off()




