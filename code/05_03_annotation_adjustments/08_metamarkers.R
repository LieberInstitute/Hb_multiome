#Getting the donor metamarkers from the current annotations


library(SingleCellExperiment)
library(MetaMarkers)
library(dplyr)
library(ggplot2)
library(qs2)
library(here)

here::here()


#Path to save any generated data
new_data_path = here('processed-data', '05_03_annotation_adjustments', '08_metamarkers')
#Path to plot directory
plot_path = here('plots','05_03_annotation_adjustments', '08_metamarkers')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

multiome_sce
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))


#Try some markers where the LHb4 inhib neurons are still just LHb4
table(multiome_sce$refined_mid_cluster)
multiome_sce$refined_noLHB4Inhib = multiome_sce$refined_mid_cluster
multiome_sce$refined_noLHB4Inhib[multiome_sce$refined_mid_cluster %in% c('Inhib_LHb_4.1','Inhib_LHb_4.2')] = 'LHb.4'
table(multiome_sce$refined_noLHB4Inhib)

#Also get markers at a broader level to test the hierarchical approach to MetaMarker annotations
class_label = rep('Non-neurons', ncol(multiome_sce))
class_label[multiome_sce$refined_mid_cluster %in% c('Inhib_LHb_4.1','Inhib_LHb_4.2','LHb.1.3.4','LHb.2.7','LHb.4',
'MHb.1','MHb.1.2', 'MHb.2', 'MHb.3')] = 'Habenula'
class_label[multiome_sce$refined_mid_cluster %in% c('Excit.Thal','Inhib.Thal')] = 'Thalamus'

multiome_sce$class_label = class_label

table(multiome_sce$class_label, multiome_sce$refined_mid_cluster)

#There's some decent donor variability in cell numbers for some clusters, so the stats might be weird, but the rankings should still be informative

all_donors = unique(multiome_sce$orig.ident)

#Get all the donor specific markers

for(i in 1:length(all_donors)){

  sce_sub = multiome_sce[, multiome_sce$orig.ident == all_donors[i]]

  markers_sub = compute_markers(assay(sce_sub, "cpm"), sce_sub$refined_mid_cluster)
  export_markers(markers_sub, paste0(new_data_path, sprintf('/markers_%s_refined_mid_markers.csv', all_donors[i])))

  markers_sub = compute_markers(assay(sce_sub, "cpm"), sce_sub$refined_noLHB4Inhib)
  export_markers(markers_sub, paste0(new_data_path, sprintf('/markers_%s_refined_noLHb4Inhib.csv', all_donors[i])))
  
  markers_sub = compute_markers(assay(sce_sub, "cpm"), sce_sub$class_label)
  export_markers(markers_sub, paste0(new_data_path, sprintf('/markers_%s_class_label.csv', all_donors[i])))


}




#Load up markers and get the metaMarkers 
multiome_refined_mid_markers = list(
    all_donors_1 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_mid_markers.csv.gz', all_donors[1]))),
    all_donors_2 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_mid_markers.csv.gz', all_donors[2]))),
    all_donors_3 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_mid_markers.csv.gz', all_donors[3]))),
    all_donors_4 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_mid_markers.csv.gz', all_donors[4]))),
    all_donors_5 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_mid_markers.csv.gz', all_donors[5]))),
    all_donors_6 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_mid_markers.csv.gz', all_donors[6]))),
    all_donors_7 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_mid_markers.csv.gz', all_donors[7]))),
    all_donors_8 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_mid_markers.csv.gz', all_donors[8]))),
    all_donors_9 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_mid_markers.csv.gz', all_donors[9]))),
    all_donors_10 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_mid_markers.csv.gz', all_donors[10])))
 
)

names(multiome_refined_mid_markers) = all_donors
multiome_refined_mid_markers

multiome_mid_metaMarkers = make_meta_markers(multiome_refined_mid_markers, detailed_stats = TRUE)

#Save the metamarkers
export_meta_markers(multiome_mid_metaMarkers, 
  paste0(new_data_path, '/multiome_refined_mid_meta_markers.csv'), 
  names(multiome_mid_metaMarkers))

multiome_mid_metaMarkers = read_meta_markers(paste0(new_data_path, '/multiome_refined_mid_meta_markers.csv.gz'))

multiome_mid_metaMarkers  %>% group_by(cell_type) %>% slice_min(rank, n = 20) %>% View()



multiome_refined_noInhib_mid_markers = list(
    all_donors_1 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_noLHb4Inhib.csv.gz', all_donors[1]))),
    all_donors_2 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_noLHb4Inhib.csv.gz', all_donors[2]))),
    all_donors_3 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_noLHb4Inhib.csv.gz', all_donors[3]))),
    all_donors_4 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_noLHb4Inhib.csv.gz', all_donors[4]))),
    all_donors_5 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_noLHb4Inhib.csv.gz', all_donors[5]))),
    all_donors_6 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_noLHb4Inhib.csv.gz', all_donors[6]))),
    all_donors_7 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_noLHb4Inhib.csv.gz', all_donors[7]))),
    all_donors_8 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_noLHb4Inhib.csv.gz', all_donors[8]))),
    all_donors_9 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_noLHb4Inhib.csv.gz', all_donors[9]))),
    all_donors_10 = read_markers(paste0(new_data_path, sprintf('/markers_%s_refined_noLHb4Inhib.csv.gz', all_donors[10])))
 
)

names(multiome_refined_noInhib_mid_markers) = all_donors
multiome_refined_noInhib_mid_markers

multiome_mid_noInhib_metaMarkers = make_meta_markers(multiome_refined_noInhib_mid_markers, detailed_stats = TRUE)

#Save the metamarkers
export_meta_markers(multiome_mid_noInhib_metaMarkers, 
  paste0(new_data_path, '/multiome_refined_mid_noInhibLHb4_meta_markers.csv'), 
  names(multiome_mid_noInhib_metaMarkers))

multiome_mid_noInhib_metaMarkers = read_meta_markers(paste0(new_data_path, '/multiome_refined_mid_noInhibLHb4_meta_markers.csv.gz'))

multiome_mid_noInhib_metaMarkers  %>% group_by(cell_type) %>% slice_min(rank, n = 20) %>% View()


#Class level metamarkers
multiome_class_markers = list(
    all_donors_1 = read_markers(paste0(new_data_path, sprintf('/markers_%s_class_label.csv.gz', all_donors[1]))),
    all_donors_2 = read_markers(paste0(new_data_path, sprintf('/markers_%s_class_label.csv.gz', all_donors[2]))),
    all_donors_3 = read_markers(paste0(new_data_path, sprintf('/markers_%s_class_label.csv.gz', all_donors[3]))),
    all_donors_4 = read_markers(paste0(new_data_path, sprintf('/markers_%s_class_label.csv.gz', all_donors[4]))),
    all_donors_5 = read_markers(paste0(new_data_path, sprintf('/markers_%s_class_label.csv.gz', all_donors[5]))),
    all_donors_6 = read_markers(paste0(new_data_path, sprintf('/markers_%s_class_label.csv.gz', all_donors[6]))),
    all_donors_7 = read_markers(paste0(new_data_path, sprintf('/markers_%s_class_label.csv.gz', all_donors[7]))),
    all_donors_8 = read_markers(paste0(new_data_path, sprintf('/markers_%s_class_label.csv.gz', all_donors[8]))),
    all_donors_9 = read_markers(paste0(new_data_path, sprintf('/markers_%s_class_label.csv.gz', all_donors[9]))),
    all_donors_10 = read_markers(paste0(new_data_path, sprintf('/markers_%s_class_label.csv.gz', all_donors[10])))
 
)

names(multiome_class_markers) = all_donors
multiome_class_markers

multiome_class_metaMarkers = make_meta_markers(multiome_class_markers, detailed_stats = TRUE)

#Save the metamarkers
export_meta_markers(multiome_class_metaMarkers, 
  paste0(new_data_path, '/multiome_class_meta_markers.csv'), 
  names(multiome_class_metaMarkers))

multiome_class_metaMarkers = read_meta_markers(paste0(new_data_path, '/multiome_class_meta_markers.csv.gz'))

multiome_class_metaMarkers  %>% group_by(cell_type) %>% slice_min(rank, n = 20) %>% View()


