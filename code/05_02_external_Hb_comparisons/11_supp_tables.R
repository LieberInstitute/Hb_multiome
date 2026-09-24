#Script to pull together all the different annotations for supplemental tables

library(SingleCellExperiment)
library(qs2)
library(dplyr)
library(here)

here::here()


#Path for new data generated
new_data_path = here('processed-data', '05_02_external_Hb_comparisons','11_supp_tables')
if (!dir.exists(new_data_path)) dir.create(new_data_path)


#Path to the Yalcinbas pilot data
#Going with the official_final_sce.RDATA
yalcinbas_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/sce_objects'
list.files(yalcinbas_path)

#Nonhuman data paths
zeb_data_path = here('processed-data', '05_02_external_Hb_comparisons','07_cross_species_hab')
mouse_data_path = here('processed-data', '05_02_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')
wallace_03_data_path = here('processed-data', '05_02_external_Hb_comparisons', '03_metaMarkers_wallace_2019')

hashikawa_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data'
hashikawa_data_path = here('processed-data','05_02_external_Hb_comparisons','04_wallace_hashikawa_mouse')


#This is the same multiome SCE as before, but with the refined annotations added in
#in refined_mid_cluster metadata
multiome_path_new = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce_new = qs_read(paste0(multiome_path_new, '/refined_annotation_multiomeHab_SCE.qs2'))

names(colData(multiome_sce_new))

#For the supplemental table, we want the fine-resolution annotations, fine-resolution annotations after minor adjustments from markers, 
# mid resolution and the consensus annotation

new_multiome_meta = as.data.frame(colData(multiome_sce_new)) |> select(cluster_ann,refined_cluster_ann, mid_cluster, refined_mid_cluster )

colnames(new_multiome_meta) = c('initial_fine_res','refined_fine_res','mid_res','consensus_annot')
#Save the metadata df as a csv file
write.csv(new_multiome_meta, paste0(new_data_path, '/multiome_annotations.csv'), row.names = TRUE)

rm(multiome_sce_new, multiome_sce)
gc()


#For all the non-multiome datasets, just need the original annotations and then the consensus annotations
#Yalcinbas data
#Loads as an object labeled 'sce'
#16437 cells
load(paste0(yalcinbas_path, '/official_final_sce.RDATA'))
yalcinbas_sce = sce
rm(sce)


names(colData(yalcinbas_sce))
#These are the original celltype annotations
table(yalcinbas_sce$final_Annotations)

meta_annot_vec = c( 'MHb_A' = 'MHb.1',
                    'MHb_B' = 'MHb.2',
                    'LHb_A' = 'LHb.7',
                    'LHb_A' = 'LHb.2',
                    'LHb_B' = 'LHb.1',
                    'LHb_B' = 'LHb.5',
                    'LHb_C' = 'LHb.3',
                    'LHb_C' = 'LHb.4',
                    'MHb_D' = 'MHb.3' ,
                    'MHb_C' = 'LHb.6',
                    'Astrocyte' = 'Astrocyte',
                    'Oligo' = 'Oligo',
                    'OPC' = 'OPC',
                    'Endo' = 'Endo',
                    'Microglia' = 'Microglia',
                    'Excit.Thal' = 'Excit.Thal',
                    'Inhib.Thal' = 'Inhib.Thal'
)

meta_annot_vec  = setNames(names(meta_annot_vec), meta_annot_vec)


yalcinbas_sce$consensus_annot = unname(meta_annot_vec[yalcinbas_sce$final_Annotations])
table(yalcinbas_sce$consensus_annot,yalcinbas_sce$final_Annotations)


yalcinbas_meta = as.data.frame(colData(yalcinbas_sce)) |> select(final_Annotations, consensus_annot)

colnames(yalcinbas_meta) = c('original_annot','consensus_annot')
#Save the metadata df as a csv file
write.csv(yalcinbas_meta, paste0(new_data_path, '/yalcinbas_annotations.csv'), row.names = TRUE)

rm(yalcinbas_sce)
gc()


#Load up the zebrafish and mouse data, filter down to the present genes
#Zebrafish data, using the summed paralog version
zeb_sce = readRDS(file = paste0(zeb_data_path, '/adult_zebrafish_summed_paralogs.rds'))
colnames(colData(zeb_sce))
#Original labels
table(zeb_sce$study_id)
#de novo clusters
table(zeb_sce$seurat_clusters)
#meta_clusters
table(zeb_sce$meta_cluster)
#Annotated metaclusters
table(zeb_sce$meta_clust_celltype_annot)

zeb_sce$grouped_annot = zeb_sce$meta_clust_celltype_annot
zeb_sce$grouped_annot[zeb_sce$meta_clust_celltype_annot %in% c('dorsolateral_left_subP')] = 'MHb_A'
zeb_sce$grouped_annot[zeb_sce$meta_clust_celltype_annot %in% c('dorsomedial_right_cholinergic','dorsomedial_right_cholinergic_2')] = 'MHb_B'
zeb_sce$grouped_annot[zeb_sce$meta_clust_celltype_annot %in% c('ventral', 'ventral_2')] = 'LHb_B'
zeb_sce$grouped_annot[zeb_sce$meta_clust_celltype_annot %in% c('inhibitory_gap43')] = 'LHb_C'


zeb_meta = as.data.frame(colData(zeb_sce)) |> select(study_id, seurat_clusters, meta_cluster, meta_clust_celltype_annot, grouped_annot)
colnames(zeb_meta) = c('original_annot','de_novo_cluster','meta_cluster','annotated_meta_cluster','consensus_annot')
zeb_meta[1:10, ]

#Save the metadata df as a csv file
write.csv(zeb_meta, paste0(new_data_path, '/zebrafish_annotations.csv'), row.names = TRUE)



#Wallace mouse data
all_mouse_sce = readRDS(paste0(mouse_data_path, '/all_donor_sce_with_denovo_clusters.rds'))

# Add the mouse metadata for wallace
current_mouse_metadata = readRDS(paste0(wallace_03_data_path, '/wallace_mouse_metaclust_celltype_annot_metadata.rds'))
colData(all_mouse_sce) = S4Vectors::DataFrame(current_mouse_metadata)

colnames(colData(all_mouse_sce))
#author annotations
table(all_mouse_sce$author_celltype)
table(all_mouse_sce$author_subHab_celltype)
#denovo clusters
table(all_mouse_sce$seurat_clusters)
#metaclusters
table(all_mouse_sce$meta_cluster)
#annotated metaclusters
table(all_mouse_sce$meta_clust_celltype_annot)

#Consensus annotation
all_mouse_sce$grouped_annot = all_mouse_sce$meta_clust_celltype_annot
all_mouse_sce$grouped_annot[all_mouse_sce$meta_clust_celltype_annot%in% c('MHb_subP')] = 'MHb_A'
all_mouse_sce$grouped_annot[all_mouse_sce$meta_clust_celltype_annot %in% c('MHb_cholinergic')] = 'MHb_B'
all_mouse_sce$grouped_annot[all_mouse_sce$meta_clust_celltype_annot %in% c('LHb_2')] = 'LHb_A'
all_mouse_sce$grouped_annot[all_mouse_sce$meta_clust_celltype_annot %in% c('LHb_1')] = 'LHb_C'

table(all_mouse_sce$grouped_annot, all_mouse_sce$meta_clust_celltype_annot)


wallace_meta = as.data.frame(colData(all_mouse_sce)) |> select(author_celltype, author_subHab_celltype,seurat_clusters, meta_cluster, meta_clust_celltype_annot, grouped_annot)
colnames(wallace_meta) = c('original_annot','original_subHab_annot','de_novo_cluster','meta_cluster','annotated_meta_cluster','consensus_annot')
wallace_meta[1:10, ]

#Save the metadata df as a csv file
write.csv(wallace_meta, paste0(new_data_path, '/wallace_mouse_annotations.csv'), row.names = TRUE)




#Hashikawa mouse data
load(paste0(hashikawa_path, '/sce_mouse_habenula.Rdata'))



hashikawa_sce_sub = sce_mouse_sub$neuron #Using just the neuron subset
rm(sce_mouse_sub)

#Add the mouse metadata
full_hashikawa_metadata = readRDS(paste0(hashikawa_data_path, '/hashikawa_mouse_neuron_metaclust_celltype_annot_metadata.rds'))
colData(hashikawa_sce_sub) = full_hashikawa_metadata 


colnames(colData(hashikawa_sce_sub))
#Original author annotations
table(hashikawa_sce_sub$celltype)
#Metacluster annotations matched to the Wallace metaclusters
table(hashikawa_sce_sub$meta_clust_celltype_annot)

#consensus annotations
hashikawa_sce_sub$grouped_annot = hashikawa_sce_sub$meta_clust_celltype_annot
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$meta_clust_celltype_annot %in% c('MHb_subP')] = 'MHb_A'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$meta_clust_celltype_annot %in% c(paste0('MHb_cholinergic_', 1:4))] = 'MHb_B'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$meta_clust_celltype_annot %in% c('LHb_2', 'LHb_1_5')] = 'LHb_A'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$meta_clust_celltype_annot %in% c('LHb_1_2', 'LHb_1_4')] = 'LHb_B'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$meta_clust_celltype_annot %in% c('LHb_1_1', 'LHb_1_3')] = 'LHb_C'

table(hashikawa_sce_sub$grouped_annot, hashikawa_sce_sub$meta_clust_celltype_annot)

hashikawa_meta = as.data.frame(colData(hashikawa_sce_sub)) |> select(celltype, meta_clust_celltype_annot, grouped_annot)
colnames(hashikawa_meta) = c('original_annot','annotated_meta_cluster','consensus_annot')
hashikawa_meta[1:10, ]

#Save the metadata df as a csv file
write.csv(hashikawa_meta, paste0(new_data_path, '/hashikawa_mouse_annotations.csv'), row.names = TRUE)


########################
#
# And now the top 50 marker gene sets
#
########################


#Human markers
deconvo_marker_path = here('processed-data','05_03_annotation_adjustments','14_deconvoBuddies_markers')
#marker_stats_MeanRatio = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_MeanRatio.rds'))
marker_stats_1vAll = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_1vAll.rds'))
#marker_stats = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_combo.rds'))

top_1vsAll_marker_df = marker_stats_1vAll %>% group_by(cellType.target) %>% filter(std.logFC.rank <= 50)
top_1vsAll_marker_df

#top_1vsAll_marker_df = marker_stats_MeanRatio %>% group_by(cellType.target) %>% filter(MeanRatio.rank <= 50)

#top_1vsAll_marker_df = marker_stats %>% group_by(cellType.target) %>% filter(MeanRatio.rank <= 50)


dup_genes = top_1vsAll_marker_df$gene[which(duplicated(top_1vsAll_marker_df$gene))]
#For each duplicate, assign it to the cell-type with the better (minimum) rank
keep_dups = top_1vsAll_marker_df %>% filter(gene %in% dup_genes) %>% group_by(gene) %>% filter(std.logFC.rank == min(std.logFC.rank))
top_1vsAll_marker_df = top_1vsAll_marker_df %>% filter(!gene %in% dup_genes)
top_1vsAll_marker_df = rbind(top_1vsAll_marker_df, keep_dups)
top_1vsAll_marker_df = top_1vsAll_marker_df %>% arrange(cellType.target)

top_1vsAll_marker_df %>% group_by(cellType.target) %>% summarise(n = n())

write.csv(top_1vsAll_marker_df, paste0(new_data_path, '/top_50_marker_genes.csv'), row.names = FALSE)





