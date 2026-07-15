#Generating the cell-type marker stats of the consensus cell-types in the hashikawa dataset

library(here)
library("DeconvoBuddies")
library(sessioninfo)
library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(Matrix)
library(edgeR)
library(MetaMarkers)
library(BSgenome.Hsapiens.UCSC.hg38)



hashikawa_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data'
hashikawa_data_path = here('processed-data','05_02_external_Hb_comparisons','04_wallace_hashikawa_mouse')

wallace_path = here('processed-data', '05_02_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')
wallace_data_path = here('processed-data', '05_02_external_Hb_comparisons', '03_metaMarkers_wallace_2019')


new_data_path = here('processed-data','17_species_diverg','01_hashikawa_markers')
dir.create(new_data_path, showWarnings = FALSE)

#Path to orthologs
path_to_orthologs = here('processed-data', '05_02_external_Hb_comparisons', 'human_mouse_zebrafish_orthologs.txt.gz')


#Wallace mouse data
wallace_sce = readRDS(paste0(wallace_path, '/all_donor_sce_with_denovo_clusters.rds'))

# Add the mouse metadata
current_mouse_metadata = readRDS(paste0(wallace_data_path, '/wallace_mouse_metaclust_celltype_annot_metadata.rds'))
colData(wallace_sce) = S4Vectors::DataFrame(current_mouse_metadata)




#Loads the sce_mouse_sub list, has two sce objects. We have the consensus annotations for just the neuronal subset object
# Add the consensus annotations, and then grab the non-neurons from the full object to piece back together

load(paste0(hashikawa_path, '/sce_mouse_habenula.Rdata'))
hashikawa_sce_all = sce_mouse_sub$all #Using just the neuron subset
hashikawa_sce_neuron = sce_mouse_sub$neuron

assay(hashikawa_sce_all, 'cpm') = MetaMarkers::convert_to_cpm(assay(hashikawa_sce_all, 'counts'))
assay(hashikawa_sce_neuron, 'cpm') = MetaMarkers::convert_to_cpm(assay(hashikawa_sce_neuron, 'counts'))

# Add in the updated mouse annotations we generated for this subset
neuron_hashikawa_metadata = readRDS(paste0(hashikawa_data_path, '/hashikawa_mouse_neuron_metaclust_celltype_annot_metadata.rds'))
#Double checking this aligns with the neuronal subset cells
table(rownames(neuron_hashikawa_metadata) == colnames(hashikawa_sce_neuron))
#Pass in the updated cell-type annotations
colData(hashikawa_sce_neuron) = neuron_hashikawa_metadata 

#Grab the non-neurons from the full dataset.
nonN_index = !grepl('Neuron',hashikawa_sce_all$celltype )
temp_nonNs_sce = hashikawa_sce_all[ , nonN_index]

#Fix up the metadata for the two objects, just swap the celltype column of the neurons for the meta_clust_celltyp_annot
hashikawa_sce_neuron$celltype = hashikawa_sce_neuron$meta_clust_celltype_annot

col_index = colnames(colData(hashikawa_sce_neuron)) != 'meta_clust_celltype_annot'
colData(hashikawa_sce_neuron) = colData(hashikawa_sce_neuron)[ , col_index]

#Fix up misaligned genes too, the neuronal subset has fewer genes
shared_genes = intersect(rownames(temp_nonNs_sce), rownames(hashikawa_sce_neuron))
temp_nonNs_sce = temp_nonNs_sce[rownames(temp_nonNs_sce) %in% shared_genes, ]
hashikawa_sce_neuron = hashikawa_sce_neuron[rownames(hashikawa_sce_neuron) %in% shared_genes, ]

#And check the order
table(rownames(temp_nonNs_sce) == rownames(hashikawa_sce_neuron))
#And combine
hashikawa_final_sce = cbind(temp_nonNs_sce, hashikawa_sce_neuron)
table(hashikawa_final_sce$celltype)


#This is the consensus cross-species mapping. Notably, we're leaving the mouse MHb_subP_cholinergic cell-type as is, as did not find a strong match in humans
#Also group the non-neuronal subtypes
hashikawa_final_sce$consensus_annot = hashikawa_final_sce$celltype
hashikawa_final_sce$consensus_annot[hashikawa_final_sce$celltype %in% c('MHb_subP')] = 'MHb.1'
hashikawa_final_sce$consensus_annot[hashikawa_final_sce$celltype %in% c(paste0('MHb_cholinergic_', 1:4))] = 'MHb.2'
hashikawa_final_sce$consensus_annot[hashikawa_final_sce$celltype %in% c('LHb_2', 'LHb_1_5')] = 'LHb.2.7'
hashikawa_final_sce$consensus_annot[hashikawa_final_sce$celltype %in% c('LHb_1_2', 'LHb_1_4')] = 'LHb.1.3.4'
hashikawa_final_sce$consensus_annot[hashikawa_final_sce$celltype %in% c('LHb_1_1', 'LHb_1_3')] = 'LHb.4'
hashikawa_final_sce$consensus_annot[grepl('Astrocyte', hashikawa_final_sce$celltype) ] = 'Astrocyte'
hashikawa_final_sce$consensus_annot[grepl('Oligo', hashikawa_final_sce$celltype) ] = 'Oligo'
hashikawa_final_sce$consensus_annot[grepl('OPC', hashikawa_final_sce$celltype) ] = 'OPC'

table(hashikawa_final_sce$consensus_annot)

#wallace annotations
wallace_sce$consensus_annot = wallace_sce$meta_clust_celltype_annot
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('MHb_subP')] = 'MHb.1'
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('MHb_cholinergic')] = 'MHb.2'
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('LHb_2')] = 'LHb.2.7'
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('LHb_1')] = 'LHb.4'
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('Astrocytes')] = 'Astrocyte'
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('Oligodendrocytes')] = 'Oligo'
wallace_sce$consensus_annot[wallace_sce$meta_clust_celltype_annot %in% c('Differentiating Oligodendrocytes')] = 'OPC'

table(wallace_sce$consensus_annot)


#Gene symbols as the rownames
rownames(hashikawa_final_sce) = rowData(hashikawa_final_sce)$Symbol
#rownames(wallace_sce) = rowData(wallace_sce)$Symbol

#Maybe clear up some memory, dont need the reduced dims for the t-stats
reducedDims(hashikawa_final_sce) <- list()
rm(temp_nonNs_sce,hashikawa_sce_neuron,hashikawa_sce_all , sce_mouse_sub  )
gc()


#And swap in the human gene name 
hu_mu_zf_ortholog_df = data.table::fread(path_to_orthologs)
table(hu_mu_zf_ortholog_df$`Mouse homology type`)

#Get the 1to1 orthologs for the mouse
hu_mu_zf_ortholog_df <- hu_mu_zf_ortholog_df %>%
  filter(
    `Mouse homology type` == 'ortholog_one2one'
  ) %>%
  filter(!duplicated(`Gene name`))
dim(hu_mu_zf_ortholog_df)

index = match(rownames(wallace_sce), hu_mu_zf_ortholog_df$`Mouse gene name`)
rownames(wallace_sce) = hu_mu_zf_ortholog_df$`Gene name`[index]

index = match(rownames(hashikawa_final_sce), hu_mu_zf_ortholog_df$`Mouse gene name`)
rownames(hashikawa_final_sce) = hu_mu_zf_ortholog_df$`Gene name`[index]

#And filter all the NAs
na_gene_index = !is.na(rownames(hashikawa_final_sce))
hashikawa_final_sce = hashikawa_final_sce[na_gene_index, ]
table(is.na(rownames(hashikawa_final_sce)))

na_gene_index = !is.na(rownames(wallace_sce))
wallace_sce = wallace_sce[na_gene_index, ]
table(is.na(rownames(wallace_sce)))


#logcounts normalization
hashikawa_final_sce = logNormCounts(hashikawa_final_sce)

#Now get the marker stats
#MeanRatio markers
marker_stats_MeanRatio <- get_mean_ratio(
    sce = hashikawa_final_sce, # sce is the SingleCellExperiment with our data
    assay_name = "logcounts", ## assay to use, we recommend logcounts [default]
    cellType_col = "consensus_annot", # column in colData with cell type info
)


#1vsAll markers
marker_stats_1vAll <- findMarkers_1vAll(
     sce = hashikawa_final_sce, # sce is the SingleCellExperiment with our data
     assay_name = "logcounts",
     cellType_col = "consensus_annot" # column in colData with cell type info
 )



## join the two marker_stats tables
marker_stats <- marker_stats_MeanRatio |>
    left_join(marker_stats_1vAll, by = join_by(gene, cellType.target))


#Save results

saveRDS(marker_stats_MeanRatio, file = paste0(new_data_path, '/marker_stats_MeanRatio.rds'))
saveRDS(marker_stats_1vAll , file = paste0(new_data_path, '/marker_stats_1vAll.rds'))
saveRDS(marker_stats, file = paste0(new_data_path, '/marker_stats_combo.rds'))


#And now generate metaMarkers too, based on the stim vs cntl labels from hashikawa, the only sample info we have
table(hashikawa_final_sce$stim)
hashikawa_cntl = hashikawa_final_sce[ , hashikawa_final_sce$stim == 'cntl' ]
hashikawa_stim = hashikawa_final_sce[ , hashikawa_final_sce$stim == 'stim' ]


#Also add in the Wallace data here too
colnames(colData(wallace_sce))
table(wallace_sce$putative_donor, wallace_sce$consensus_annot)
wallace_d1 = wallace_sce[ , wallace_sce$putative_donor == '160822' ]
wallace_d2 = wallace_sce[ , wallace_sce$putative_donor == '161102' ]
wallace_d3 = wallace_sce[ , wallace_sce$putative_donor == '161103' ]
wallace_d4 = wallace_sce[ , wallace_sce$putative_donor == '161105' ]

mouse_marker_list = list(
  cntl = compute_markers(assay(hashikawa_cntl, 'cpm'), hashikawa_cntl$consensus_annot),
  stim = compute_markers(assay(hashikawa_stim, 'cpm'), hashikawa_stim$consensus_annot),
  wallace_d1 = compute_markers(assay(wallace_d1, 'cpm'), wallace_d1$consensus_annot),
  wallace_d2 = compute_markers(assay(wallace_d2, 'cpm'), wallace_d2$consensus_annot),
  wallace_d3 = compute_markers(assay(wallace_d3, 'cpm'), wallace_d3$consensus_annot),
  wallace_d4 = compute_markers(assay(wallace_d4, 'cpm'), wallace_d4$consensus_annot)
  )


mouse_metaMarkers = make_meta_markers(mouse_marker_list , detailed_stats = TRUE)

#Save the metamarkers
export_meta_markers(mouse_metaMarkers, 
  paste0(new_data_path, '/mouse_meta_markers.csv'), 
  names(mouse_metaMarkers))


session_info()