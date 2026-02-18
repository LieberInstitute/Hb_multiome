#Script for checking out the MetaMarkers across the zebrafish samples
#And initial exploration of the same markers used to annotate the mouse clusters



library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
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
new_data_path = here('processed-data', '98_external_Hb_comparisons', '05_metaMarkers_pandey_2018')
#Path to already generated data
prev_data_path = here('processed-data', '98_external_Hb_comparisons', '04_initial_qc_pandey_2018')
#Path to plot directory
plot_path = here('plots', '98_external_Hb_comparisons', '05_metaMarkers_pandey_2018')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)



#Load up the full SCE object, contains the metacluster annotations
all_donor_sce = readRDS(paste0(prev_data_path, '/all_donor_sce_with_denovo_clusters.rds'))

#Generate MetaMarkers from the metaclusters
#Split into individual donor SCE objects and compute intial markers for the metaclusters
table(all_donor_sce$study_id)

sce_zebAdult1 <- all_donor_sce[, all_donor_sce$study_id == "Adult1"]
sce_zebAdult2 <- all_donor_sce[, all_donor_sce$study_id == "Adult2"]
sce_zebLarva <- all_donor_sce[, all_donor_sce$study_id == "Larva"]


markers_adult1 = compute_markers(assay(sce_zebAdult1, "cpm"), sce_zebAdult1$meta_cluster)
markers_adult2 = compute_markers(assay(sce_zebAdult2, "cpm"), sce_zebAdult2$meta_cluster)
markers_larva = compute_markers(assay(sce_zebLarva, "cpm"), sce_zebLarva$meta_cluster)

head(markers_adult1)


#Save markers
export_markers(markers_adult1, paste0(new_data_path, '/markers_adult1_meta_clusters_markers.csv'))
export_markers(markers_adult2, paste0(new_data_path, '/markers_adult2_meta_clusters_markers.csv'))
export_markers(markers_larva, paste0(new_data_path, '/markers_larva_meta_clusters_markers.csv'))


#Load up markers and get the metaMarkers 
pandey_zebrafish_hab_markers = list(
    markers_adult1 = read_markers(paste0(new_data_path, '/markers_adult1_meta_clusters_markers.csv.gz')),
    markers_adult2 = read_markers(paste0(new_data_path, '/markers_adult2_meta_clusters_markers.csv.gz')),
    markers_larva = read_markers(paste0(new_data_path, '/markers_larva_meta_clusters_markers.csv.gz'))
    
)


pandey_zebrafish_hab_metaM = make_meta_markers(pandey_zebrafish_hab_markers, detailed_stats = TRUE)

#Save the metamarkers
export_meta_markers(pandey_zebrafish_hab_metaM, 
  paste0(new_data_path, '/zebrafish_hab_meta_markers.csv'), 
  names(pandey_zebrafish_hab_metaM))

pandey_zebrafish_hab_metaM = read_meta_markers(paste0(new_data_path, '/zebrafish_hab_meta_markers.csv.gz'))

pandey_zebrafish_hab_metaM %>% group_by(cell_type) %>% slice_min(rank, n = 20) %>% View()















