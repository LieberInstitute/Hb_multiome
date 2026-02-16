#Script for getting MetaMarkers from the metaClusters across individual mouse donors from the Wallace 2019 dataset
#Will likely be updating metaCluster annotations, just setting up the code to check out the markers for those clusters


library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
library(dplyr)
library(ggplot2)




#Load up the full SCE object, contains the metacluster annotations
all_donor_sce = readRDS('processed-data/98_external_Hb_comparisons/02_qc_and_clust_wallace_2019/all_donor_sce_with_denovo_clusters.rds')


#Generate MetaMarkers from the metaclusters
#Split into individual donor SCE objects and compute intial markers for the metaclusters
table(all_donor_sce$putative_donor)

sce_160822 <- all_donor_sce[, all_donor_sce$putative_donor == "160822"]
sce_161102 <- all_donor_sce[, all_donor_sce$putative_donor == "161102"]
sce_161103 <- all_donor_sce[, all_donor_sce$putative_donor == "161103"]
sce_161105 <- all_donor_sce[, all_donor_sce$putative_donor == "161105"]

markers_160822 = compute_markers(assay(sce_160822, "cpm"), sce_160822$meta_cluster)
markers_161102 = compute_markers(assay(sce_161102, "cpm"), sce_161102$meta_cluster)
markers_161103 = compute_markers(assay(sce_161103, "cpm"), sce_161103$meta_cluster)
markers_161105 = compute_markers(assay(sce_161105, "cpm"), sce_161105$meta_cluster)

head(markers_160822)

#Save markers
export_markers(markers_160822, 'processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_160822_meta_clusters_markers.csv')
export_markers(markers_161103, 'processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_161103_meta_clusters_markers.csv')
export_markers(markers_161102, 'processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_161102_meta_clusters_markers.csv')
export_markers(markers_161105, 'processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_161105_meta_clusters_markers.csv')


#Load up markers and get the metaMarkers 
wallace_mouse_hab_markers = list(
    donor_160822 = read_markers("processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_160822_meta_clusters_markers.csv.gz"),
    donor_161102 = read_markers("processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_161102_meta_clusters_markers.csv.gz"),
    donor_161103 = read_markers("processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_161103_meta_clusters_markers.csv.gz"),
    donor_161105 = read_markers("processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/markers_161105_meta_clusters_markers.csv.gz")
    
)


wallace_mouse_hab_metaM = make_meta_markers(wallace_mouse_hab_markers, detailed_stats = TRUE)

#Save the metamarkers
export_meta_markers(wallace_mouse_hab_metaM, 
  "processed-data/98_external_Hb_comparisons/03_metaMarkers_wallace_2019/mouse_hab_meta_markers.csv", 
  names(wallace_mouse_hab_metaM))

wallace_mouse_hab_metaM %>% group_by(cell_type) %>% slice_min(rank, n = 10) %>% View()


#Load up the individual seurat objects, add the metaCluster annotations, and check out what it looks like in the UMAPs
indv_mouse_seurats = readRDS(file = 'processed-data/98_external_Hb_comparisons/02_qc_and_clust_wallace_2019/individual_donor_seurat_objects_list.rds')

hab_160822_seurat = indv_mouse_seurats$hab_160822
hab_161102_seurat = indv_mouse_seurats$hab_161102
hab_161103_seurat = indv_mouse_seurats$hab_161103
hab_161105_seurat = indv_mouse_seurats$hab_161105


full_cell_names_160822 = paste('hab_160822', hab_160822_seurat$seurat_clusters, sep = '|')
hab_160822_seurat$meta_cluster = unname(all_donor_sce$meta_cluster[full_cell_names_160822])

full_cell_names_161102 = paste('hab_161102', hab_161102_seurat$seurat_clusters, sep = '|')
hab_161102_seurat$meta_cluster = unname(all_donor_sce$meta_cluster[full_cell_names_161102])

full_cell_names_161103 = paste('hab_161103', hab_161103_seurat$seurat_clusters, sep = '|')
hab_161103_seurat$meta_cluster = unname(all_donor_sce$meta_cluster[full_cell_names_161103])

full_cell_names_161105 = paste('hab_161105', hab_161105_seurat$seurat_clusters, sep = '|')
hab_161105_seurat$meta_cluster = unname(all_donor_sce$meta_cluster[full_cell_names_161105])



DimPlot(hab_160822_seurat, group.by = 'meta_cluster', label = TRUE) + ggtitle('Mouse: 160822')
DimPlot(hab_161102_seurat, group.by = 'meta_cluster', label = TRUE) + ggtitle('Mouse: 161102')
DimPlot(hab_161103_seurat, group.by = 'meta_cluster', label = TRUE) + ggtitle('Mouse: 161103')
DimPlot(hab_161105_seurat, group.by = 'meta_cluster', label = TRUE) + ggtitle('Mouse: 161105')



















