#Script for initial basic QC and attempting default Seurat clustering of the Wallace 2019 mouse data.
#Also add in the author provided cluster annotations

library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(MetaNeighbor)
library(ComplexHeatmap)
library(here)

here::here()

#Path to the Wallace 2019 data
wallace_path = here('processed-data', '05_02_external_Hb_comparisons', 'Wallace_etal_2019_habenula_scseq')
list.files(wallace_path)
#Path to save any generated data
new_data_path = here('processed-data', '05_02_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')
#Path to plot directory
plot_path = here('plots', '05_02_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


#Rdata objects containing seurat objects with cluster annotations from the original paper.
#As in Figure 1B
load(paste0(wallace_path, '/Habenula_Seurat_meta.RData'))
View(meta)
table(meta$CellClassNames_filtered)
dim(meta)


#Rdata objects with the Lateral and Medial habenula subclustering
#From https://github.com/mwall2017/habenula_indrops
#Cell type annotations are in the tree.ident column
#Medial Hab: 1= Ventral 2/3, 2= Ventrolateral, 3= Lateral, 4= Dorsal, 5= Superior

load(paste0(wallace_path, '/mhb_Seurat_meta.RData'))
View(meta_mhb)

dim(meta_mhb)
table(rownames(meta_mhb) %in% rownames(meta))

mhb_labels = c('Mhb_Ventral 2/3' = '1', 'Mhb_Ventrolateral' = '2', 'Mhb_Lateral' = '3', 'Mhb_Dorsal' = '4', 'Mhb_Superior' = '5')
meta_mhb$celltype_annot = names(mhb_labels[meta_mhb$tree.ident])
table(meta_mhb$celltype_annot)

#Lateral Hab: 1= Oval/Medial, 2= Marginal, 3= Lateral, 4= Hbx 
load(paste0(wallace_path, '/lhb_Seurat_meta.RData'))
View(meta_lhb)

dim(meta_lhb)
table(rownames(meta_lhb) %in% rownames(meta))

lhb_labels = c('Lhb_Oval/Medial' = '1', 'Lhb_Marginal' = '2', 'Lhb_Lateral' = '3', 'Lhb_Hox' = '4')
meta_lhb$celltype_annot = names(lhb_labels[meta_lhb$tree.ident])
table(meta_lhb$celltype_annot)

#Add to the full metadata dataframe
#Initial check separately to see what the habenula subclusters map to the original annoations
test_index = match(rownames(meta_lhb), rownames(meta))
table(meta$CellClassNames_filtered[test_index])
test_index = match(rownames(meta_mhb), rownames(meta))
table(meta$CellClassNames_filtered[test_index])
#they're not a perfect match between medial habenula in the original and then medial habenula in the subclusters
#but pretty close

mhb_index = match(rownames(meta_mhb), rownames(meta))
lhb_index = match(rownames(meta_lhb), rownames(meta))

meta$hab_subcluster = 'Not habenula'
meta$hab_subcluster[mhb_index] = meta_mhb$celltype_annot
meta$hab_subcluster[lhb_index] = meta_lhb$celltype_annot

#Confusion matrix to visualize annotation mappings
cell_num_confs_mat = as.matrix(table(meta$CellClassNames_filtered, meta$hab_subcluster))
subclust_cell_nums = colSums(cell_num_confs_mat)
cell_num_confs_mat = sweep(cell_num_confs_mat, 2, subclust_cell_nums, "/")

col_fun = circlize::colorRamp2(c(0, 1), c("white", "red"))
Heatmap(cell_num_confs_mat, name = 'Proportion of cells', col = col_fun, column_title = 'Habenula subcluster mapping' ,
cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE, show_column_names = TRUE )



#Load the data
#Cell count matches the reported 7506 final cells used in the paper
hab_batch1_data = readRDS(paste0(wallace_path, '/hab_batch1.rds'))
hab_batch1_data

#Made with Seurat v2.3.4
hab_batch1_data@version

#Update the object to v3
hab_batch1_data <- UpdateSeuratObject(hab_batch1_data)


#Check that the cluster metadata and the cell barcodes in the seurat object match up
#They do and they're in the correct order already, neat
table(rownames(hab_batch1_data@meta.data) == rownames(meta))

#Add the author cluster annotations
hab_batch1_data$author_celltype = meta$CellClassNames_filtered
hab_batch1_data$author_subHab_celltype = meta$hab_subcluster

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
DimPlot(hab_batch1_data, group.by = "author_celltype", label= TRUE)
DimPlot(hab_batch1_data, group.by = "author_subHab_celltype", label= TRUE)
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
  donor_list[[i]] <- FindClusters(donor_list[[i]], resolution = .8)
  
  # UMAP
  donor_list[[i]] <- RunUMAP(donor_list[[i]], dims = 1:20)
  
  cat("Finished processing", names(donor_list)[i], "\n")
}

#Check out the UMAPs for each donor
for (i in seq_along(donor_list)) {
  p <- DimPlot(donor_list[[i]], reduction = "umap", label = TRUE) + ggtitle(names(donor_list)[i])
  print(p)
  p <- DimPlot(donor_list[[i]], reduction = "umap",group.by = 'author_celltype' , label = TRUE) + ggtitle(names(donor_list)[i])
  print(p)
  p <- DimPlot(donor_list[[i]], reduction = "umap",group.by = 'author_subHab_celltype' , label = TRUE) + ggtitle(names(donor_list)[i])
  print(p)
}


#Save the donor-specific Seurat objects for future use
saveRDS(donor_list, file = paste0(new_data_path, '/individual_donor_seurat_objects_list.rds'))



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


#Confusion matrix between metaclusters and the author annotations
all_celltype_conf_mat = as.matrix(table(all_donor_sce$meta_cluster, all_donor_sce$author_celltype))
all_celltype_sum_vec = colSums(all_celltype_conf_mat)
all_celltype_conf_mat  = sweep(all_celltype_conf_mat , 2, all_celltype_sum_vec, "/")

col_fun = circlize::colorRamp2(c(0, 1), c("white", "red"))
Heatmap(all_celltype_conf_mat, name = 'Proportion of cells', col = col_fun, column_title = 'Metacluster vs author annotation' ,
cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE )

#Subclust habenula annotations
all_celltype_conf_mat = as.matrix(table(all_donor_sce$meta_cluster, all_donor_sce$author_subHab_celltype))
all_celltype_sum_vec = colSums(all_celltype_conf_mat)
all_celltype_conf_mat  = sweep(all_celltype_conf_mat , 2, all_celltype_sum_vec, "/")

col_fun = circlize::colorRamp2(c(0, 1), c("white", "red"))
Heatmap(all_celltype_conf_mat, name = 'Proportion of cells', col = col_fun, column_title = 'Metacluster vs author Hab-subclust annots' ,
cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE )



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
  title = "MetaNeighbor AUROCs for Wallace 2019 donors")

dev.off()


pdf(paste0(plot_path, '/best_vs_next_MN_aurocs_heatmap.pdf'), width = 10, height = 8)
#Plot best_vs_next AUROC heatmap
plotHeatmap(MN_best_aurocs, 
  show_dendro = TRUE, 
  show_labels = TRUE, 
  cex = .5,
  title = "MetaNeighbor best_vs_next AUROCs for Wallace 2019 donors") 
dev.off()




