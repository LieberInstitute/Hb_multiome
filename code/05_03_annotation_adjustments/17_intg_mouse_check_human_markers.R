#For a focused cross-species comparison
#Integrate the two mouse datasets together, plot the grouped cross-species clusters we've defined
#And then visualize the top X human markers for cell-type within the mouse data.

#A more direct comparisons than just a marker gene bubble plot across species


library(SingleCellExperiment)
library(harmony)
library(MetaMarkers)
library(dplyr)
library(tidyr)
library(ggplot2)
library(qs2)
library(scran)
library(scater)
library(here)

here::here()


plot_path = here('plots','05_03_annotation_adjustments', '17_intg_mouse_check_human_markers')
new_data_path = here('processed-data','05_03_annotation_adjustments', '17_intg_mouse_check_human_markers')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


#Human markers
deconvo_marker_path = here('processed-data','05_03_annotation_adjustments','14_deconvoBuddies_markers')
marker_stats_MeanRatio = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_MeanRatio.rds'))
marker_stats_1vAll = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_1vAll.rds'))
marker_stats = readRDS(file = paste0(deconvo_marker_path, '/marker_stats_combo.rds'))


#First load up and integrate the mouse datasets
mouse_data_path = here('processed-data', '05_02_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')
wallace_03_data_path = here('processed-data', '05_02_external_Hb_comparisons', '03_metaMarkers_wallace_2019')

hashikawa_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data'
hashikawa_data_path = here('processed-data','05_02_external_Hb_comparisons','04_wallace_hashikawa_mouse')


#Wallace mouse data
all_mouse_sce = readRDS(paste0(mouse_data_path, '/all_donor_sce_with_denovo_clusters.rds'))
all_mouse_sce

#Hashikawa mouse data
load(paste0(hashikawa_path, '/sce_mouse_habenula.Rdata'))
hashikawa_sce_sub = sce_mouse_sub$all

#Gene symbols as the rownames
rownames(hashikawa_sce_sub) = rowData(hashikawa_sce_sub)$Symbol

#Add CPM
assay(hashikawa_sce_sub, "cpm") = MetaMarkers::convert_to_cpm(assay(hashikawa_sce_sub, "counts"))


#Intersect to shared genes
keep_genes = intersect(rownames(all_mouse_sce), rownames(hashikawa_sce_sub))
length(keep_genes)

all_mouse_sce = all_mouse_sce[rownames(all_mouse_sce) %in% keep_genes, ]
hashikawa_sce_sub = hashikawa_sce_sub[rownames(hashikawa_sce_sub) %in% keep_genes, ]

# After subsetting to shared genes, reorder to match
hashikawa_sce_sub = hashikawa_sce_sub[rownames(all_mouse_sce), ]
#Check
table(rownames(all_mouse_sce) == rownames(hashikawa_sce_sub))


#Add in the current metacluster annotations

# Add the mouse metadata for wallace
current_mouse_metadata = readRDS(paste0(wallace_03_data_path, '/wallace_mouse_metaclust_celltype_annot_metadata.rds'))
#Check
table(rownames(colData(all_mouse_sce)) == rownames(current_mouse_metadata))
all_mouse_sce$meta_clust_celltype_annot = current_mouse_metadata$meta_clust_celltype_annot


hashikawa_metadata = readRDS(paste0(hashikawa_data_path, '/hashikawa_mouse_metaclust_celltype_annot_metadata.rds'))
table(rownames(colData(hashikawa_sce_sub)) == rownames(hashikawa_metadata))
hashikawa_sce_sub$meta_clust_celltype_annot = hashikawa_metadata$meta_clust_celltype_annot


neuronal_subset_metadata = readRDS(paste0(hashikawa_data_path, '/hashikawa_mouse_neuron_metaclust_celltype_annot_metadata.rds'))

#Pass in the neuronal subset annotations
index = match(rownames(neuronal_subset_metadata), colnames(hashikawa_sce_sub) )
hashikawa_sce_sub$meta_clust_celltype_annot[index] = neuronal_subset_metadata$meta_clust_celltype_annot

#Give a generic neuron label for the cells that were not included in the neuronal subset
table(hashikawa_sce_sub$meta_clust_celltype_annot)
index = hashikawa_sce_sub$meta_clust_celltype_annot %in% c('LHb_1.1','LHb_1.2','LHb_1.3','MHb.1','MHb.2','MHb.3', 'MHb_cholinergic.1','MHb_cholinergic.2','MHb_cholinergic.3')
hashikawa_sce_sub$meta_clust_celltype_annot[index] = 'Neuron'




#And the annotations reflecting the cross-species mapping
all_mouse_sce$grouped_annot = all_mouse_sce$meta_clust_celltype_annot
all_mouse_sce$grouped_annot[all_mouse_sce$meta_clust_celltype_annot %in% c('MHb_subP')] = 'MHb.1'
all_mouse_sce$grouped_annot[all_mouse_sce$meta_clust_celltype_annot %in% c('MHb_subP_cholinergic')] = 'MHb.1.2'
all_mouse_sce$grouped_annot[all_mouse_sce$meta_clust_celltype_annot %in% c('MHb_cholinergic')] = 'MHb.2'
all_mouse_sce$grouped_annot[all_mouse_sce$meta_clust_celltype_annot %in% c('LHb_2')] = 'LHb.2.7'
all_mouse_sce$grouped_annot[all_mouse_sce$meta_clust_celltype_annot %in% c('LHb_1')] = 'LHb.1.3.4'
all_mouse_sce$grouped_annot[all_mouse_sce$meta_clust_celltype_annot %in% c('Astrocytes', 'Differentiating Oligodendrocytes',
'Endothelial','Fibroblasts','Macrophages','Microglia', 'Oligodendrocytes', 'Pericytes','Polydendrocytes')] = 'Non-neurons'

hashikawa_sce_sub$grouped_annot = hashikawa_sce_sub$meta_clust_celltype_annot
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$meta_clust_celltype_annot %in% c('MHb_subP')] = 'MHb.1'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$meta_clust_celltype_annot %in% c('MHb_subP_cholinergic')] = 'MHb.1.2'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$meta_clust_celltype_annot %in% c(paste0('MHb_cholinergic_', 1:4))] = 'MHb.2'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$meta_clust_celltype_annot %in% c('LHb_2', 'LHb_1_5')] = 'LHb.2.7'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$meta_clust_celltype_annot %in% c('LHb_1_2', 'LHb_1_4')] = 'LHb.1.3.4'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$meta_clust_celltype_annot %in% c('LHb_1_1', 'LHb_1_3')] = 'LHb.4'
hashikawa_sce_sub$grouped_annot[hashikawa_sce_sub$meta_clust_celltype_annot %in% c('Astrocyte1', 'Astrocyte2', 'Endothelial', 'Epen', 'Microglia', 
'Mural', 'Oligo1', 'Oligo2', 'Oligo3', 'OPC1', 'OPC2', 'OPC3')] = 'Non-neurons'
table(all_mouse_sce$grouped_annot)



#Variable to run correction on will be putative_donor for wallace, and stim for hashikawa

all_mouse_sce$batch_variable = all_mouse_sce$putative_donor
hashikawa_sce_sub$batch_variable = hashikawa_sce_sub$stim

all_mouse_sce$batch_variable_2 = 'Wallace'
hashikawa_sce_sub$batch_variable_2 = 'Hashikawa'


all_mouse_sce$original_annotations = all_mouse_sce$author_subHab_celltype
hashikawa_sce_sub$original_annotations = hashikawa_sce_sub$celltype


assay(all_mouse_sce, "scaledata") = NULL

keep_index = colnames(colData(all_mouse_sce)) %in% c('batch_variable','batch_variable_2', 'meta_clust_celltype_annot', 'original_annotations', 'grouped_annot')
colData(all_mouse_sce) = colData(all_mouse_sce)[, keep_index]

keep_index = colnames(colData(hashikawa_sce_sub)) %in% c('batch_variable', 'batch_variable_2', 'meta_clust_celltype_annot', 'original_annotations', 'grouped_annot')
colData(hashikawa_sce_sub) = colData(hashikawa_sce_sub)[, keep_index]


#Do each dataset separately
all_mouse_sce = all_mouse_sce[, !all_mouse_sce$grouped_annot %in% c('outliers', 'Neuron', 'Non-neurons')]
all_mouse_sce = logNormCounts(all_mouse_sce,name = "logcounts")

dec.all_mouse_sce <- modelGeneVar(all_mouse_sce, assay.type = "logcounts")
hvgs <- getTopHVGs(dec.all_mouse_sce, n = 2000)

all_mouse_sce <- runPCA(all_mouse_sce,
    subset_row = hvgs,
    ncomponents = 30,
    name = "PCA"
)

message("Running TSNE - ", Sys.time())
all_mouse_sce <- runTSNE(all_mouse_sce, dimred = "PCA", name = "TSNE.PCA")
colnames(reducedDim(all_mouse_sce, "TSNE.PCA")) <- c("TSNE1", "TSNE2")

message("Running UMAP - ", Sys.time())
all_mouse_sce <- runUMAP(all_mouse_sce, dimred = "PCA", name = "UMAP.PCA")
colnames(reducedDim(all_mouse_sce, "UMAP.PCA")) <- c("UMAP1", "UMAP2")


#Hashikawa dataset
hashikawa_sce_sub = hashikawa_sce_sub[, !hashikawa_sce_sub$grouped_annot %in% c('outliers', 'Neuron', 'Non-neurons')]
hashikawa_sce_sub = logNormCounts(hashikawa_sce_sub,name = "logcounts")

dec.hashikawa_sce_sub <- modelGeneVar(hashikawa_sce_sub, assay.type = "logcounts")
hvgs <- getTopHVGs(dec.hashikawa_sce_sub, n = 2000)

hashikawa_sce_sub <- runPCA(hashikawa_sce_sub,
    subset_row = hvgs,
    ncomponents = 30,
    name = "PCA"
)

message("Running TSNE - ", Sys.time())
hashikawa_sce_sub <- runTSNE(hashikawa_sce_sub, dimred = "PCA", name = "TSNE.PCA")
colnames(reducedDim(hashikawa_sce_sub, "TSNE.PCA")) <- c("TSNE1", "TSNE2")

message("Running UMAP - ", Sys.time())
hashikawa_sce_sub <- runUMAP(hashikawa_sce_sub, dimred = "PCA", name = "UMAP.PCA")
colnames(reducedDim(hashikawa_sce_sub, "UMAP.PCA")) <- c("UMAP1", "UMAP2")




#both_mouse_sce = cbind(all_mouse_sce, hashikawa_sce_sub)

#Filter out the outlier and Neuron cells, and the non-neurons
#both_mouse_sce = both_mouse_sce[, !both_mouse_sce$grouped_annot %in% c('outliers', 'Neuron', 'Non-neurons')]


#both_mouse_sce = logNormCounts(both_mouse_sce,name = "logcounts")

#And standard HVG
#dec.both_mouse_sce <- modelGeneVar(both_mouse_sce, assay.type = "logcounts")
#hvgs <- getTopHVGs(dec.both_mouse_sce, n = 2000)


#both_mouse_sce <- runPCA(both_mouse_sce,
#    subset_row = hvgs,
#    ncomponents = 30,
#    name = "PCA"
#)




# ## Run Harmony
# message("Running Harmony - ", Sys.time())
# both_mouse_sce <- RunHarmony(both_mouse_sce, group.by.vars = c('batch_variable_2'), verbose = TRUE)


# #Visualize before and after correction
# message("Running TSNE - ", Sys.time())
# both_mouse_sce <- runTSNE(both_mouse_sce, dimred = "HARMONY", name = "TSNE.HARMONY")
# colnames(reducedDim(both_mouse_sce, "TSNE.HARMONY")) <- c("TSNE1", "TSNE2")

# message("Running UMAP - ", Sys.time())
# both_mouse_sce <- runUMAP(both_mouse_sce, dimred = "HARMONY", name = "UMAP.HARMONY")
# colnames(reducedDim(both_mouse_sce, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")


# message("Running TSNE - ", Sys.time())
# both_mouse_sce <- runTSNE(both_mouse_sce, dimred = "PCA", name = "TSNE.PCA")
# colnames(reducedDim(both_mouse_sce, "TSNE.PCA")) <- c("TSNE1", "TSNE2")

# message("Running UMAP - ", Sys.time())
# both_mouse_sce <- runUMAP(both_mouse_sce, dimred = "PCA", name = "UMAP.PCA")
# colnames(reducedDim(both_mouse_sce, "UMAP.PCA")) <- c("UMAP1", "UMAP2")



plotReducedDim(
  hashikawa_sce_sub,
  dimred = "UMAP.PCA", colour_by = "batch_variable"
)

plotReducedDim(
  hashikawa_sce_sub,
  dimred = "TSNE.PCA", colour_by = "grouped_annot"
)


plotReducedDim(
  all_mouse_sce,
  dimred = "UMAP.PCA", colour_by = "batch_variable"
)

plotReducedDim(
  all_mouse_sce,
  dimred = "TSNE.PCA", colour_by = "grouped_annot"
)


#And now visualize the top human markers in this integrated data.


#Use the MetaMarkers framework to get marker set enrichments per cell-type per cell
#Then visualize those on the umap


# Add in the human gene names to the mouse data

#Path to orthologs
path_to_orthologs = here('processed-data', '05_02_external_Hb_comparisons', 'human_mouse_zebrafish_orthologs.txt.gz')

hu_mu_zf_ortholog_df = data.table::fread(path_to_orthologs)

#Get the 1to1 orthologs for the mouse
hu_mu_zf_ortholog_df <- hu_mu_zf_ortholog_df %>%
  filter(
    `Mouse homology type` == 'ortholog_one2one'
  ) %>%
  filter(!duplicated(`Gene name`))
dim(hu_mu_zf_ortholog_df)


index = match(rownames(all_mouse_sce), hu_mu_zf_ortholog_df$`Mouse gene name`)
rowData(all_mouse_sce)$human_gene_name = hu_mu_zf_ortholog_df$`Gene name`[index]

ortho_present_mouse_sce = all_mouse_sce[!is.na(rowData(all_mouse_sce)$human_gene_name) , ]
rownames(ortho_present_mouse_sce) = rowData(ortho_present_mouse_sce)$human_gene_name


index = match(rownames(hashikawa_sce_sub), hu_mu_zf_ortholog_df$`Mouse gene name`)
rowData(hashikawa_sce_sub)$human_gene_name = hu_mu_zf_ortholog_df$`Gene name`[index]

ortho_present_hashikawa_sce_sub = hashikawa_sce_sub[!is.na(rowData(hashikawa_sce_sub)$human_gene_name) , ]
rownames(ortho_present_hashikawa_sce_sub) = rowData(ortho_present_hashikawa_sce_sub)$human_gene_name



original_metadata_wallace = colData(ortho_present_mouse_sce)
original_metadata_hashikawa = colData(ortho_present_hashikawa_sce_sub)




top_1vsAll_marker_df = marker_stats_1vAll %>% group_by(cellType.target) %>% filter(std.logFC.rank <= 50)
top_1vsAll_marker_df

top_1vsAll_marker_df = marker_stats_MeanRatio %>% group_by(cellType.target) %>% filter(MeanRatio.rank <= 50)

top_1vsAll_marker_df = marker_stats %>% group_by(cellType.target) %>% filter(MeanRatio.rank <= 50)


dup_genes = top_1vsAll_marker_df$gene[which(duplicated(top_1vsAll_marker_df$gene))]
#For each duplicate, assign it to the cell-type with the better (minimum) rank
keep_dups = top_1vsAll_marker_df %>% filter(gene %in% dup_genes) %>% group_by(gene) %>% filter(std.logFC.rank == min(std.logFC.rank))
top_1vsAll_marker_df = top_1vsAll_marker_df %>% filter(!gene %in% dup_genes)
top_1vsAll_marker_df = rbind(top_1vsAll_marker_df, keep_dups)
top_1vsAll_marker_df = top_1vsAll_marker_df %>% arrange(cellType.target)

top_1vsAll_marker_df %>% group_by(cellType.target) %>% summarise(n = n())





#Markers present in the Wallace data
top_current_markers_all_mouse = top_1vsAll_marker_df %>% filter(cellType.target %in% c('MHb.1', 'MHb.1.2', 'MHb.2', 'MHb.3', 'LHb.1.3.4', 'LHb.2.7','LHb.4') ) %>%
  select(gene, cellType.target) %>% filter(gene %in% rowData(ortho_present_mouse_sce)$human_gene_name)

colnames(top_current_markers_all_mouse) = c('gene', 'cell_type')
top_current_markers_all_mouse$group = 'All'

#Markers present in the Hashikawa dataset
top_current_markers_hashikawa = top_1vsAll_marker_df %>% filter(cellType.target %in% c('MHb.1', 'MHb.1.2', 'MHb.2', 'MHb.3', 'LHb.1.3.4', 'LHb.2.7','LHb.4') ) %>%
  select(gene, cellType.target) %>% filter(gene %in% rowData(ortho_present_hashikawa_sce_sub)$human_gene_name)

colnames(top_current_markers_hashikawa) = c('gene', 'cell_type')
top_current_markers_hashikawa$group = 'All'



#Marker enrichment Wallace
ct_scores = score_cells(log1p(cpm(ortho_present_mouse_sce)), top_current_markers_all_mouse)
ct_enrichment = compute_marker_enrichment(ct_scores)

scale(t(ct_enrichment))

# Add in the enrichments to the metadata

colData(ortho_present_mouse_sce) = cbind(original_metadata_wallace, scale(t(ct_enrichment)))


#Marker enrichment Hashikawa
ct_scores = score_cells(log1p(cpm(ortho_present_hashikawa_sce_sub)), top_current_markers_hashikawa)
ct_enrichment = compute_marker_enrichment(ct_scores)

# Add in the enrichments to the metadata
colData(ortho_present_hashikawa_sce_sub) = cbind(original_metadata_hashikawa,scale(t(ct_enrichment)))



plotReducedDim(
  ortho_present_mouse_sce,
  dimred = "TSNE.PCA", colour_by = "grouped_annot"
)

plotReducedDim(
  ortho_present_mouse_sce,
  dimred = "TSNE.PCA", colour_by = "All|MHb.1"
) + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(NA, 4), oob = scales::squish) +
  ggtitle('Human MHb.1 markers in Wallace')

plotReducedDim(
  ortho_present_mouse_sce,
  dimred = "TSNE.PCA", colour_by = "All|MHb.1.2"
) + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(NA, 4), oob = scales::squish) +
  ggtitle('Human MHb.1.2 markers in Wallace')

plotReducedDim(
  ortho_present_mouse_sce,
  dimred = "TSNE.PCA", colour_by = "All|MHb.2"
) + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(NA, 4), oob = scales::squish) +
  ggtitle('Human MHb.2 markers in Wallace')

plotReducedDim(
  ortho_present_mouse_sce,
  dimred = "TSNE.PCA", colour_by = "All|LHb.1.3.4"
) + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(NA, 4), oob = scales::squish) +
  ggtitle('Human LHb.1.3.4 markers in Wallace')

plotReducedDim(
  ortho_present_mouse_sce,
  dimred = "TSNE.PCA", colour_by = "All|LHb.2.7"
) + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(NA, 4), oob = scales::squish) +
  ggtitle('Human LHb.2.7 markers in Wallace')

plotReducedDim(
  ortho_present_mouse_sce,
  dimred = "TSNE.PCA", colour_by = "All|LHb.4"
) + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(NA, 4), oob = scales::squish) +
  ggtitle('Human LHb.4 markers in Wallace')






#Hashikawa
plotReducedDim(
  ortho_present_hashikawa_sce_sub,
  dimred = "TSNE.PCA", colour_by = "grouped_annot"
)

plotReducedDim(
  ortho_present_hashikawa_sce_sub,
  dimred = "TSNE.PCA", colour_by = "All|MHb.1"
) + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(NA, 4), oob = scales::squish) +
  ggtitle('Human MHb.1 markers in Hashikawa')

plotReducedDim(
  ortho_present_hashikawa_sce_sub,
  dimred = "TSNE.PCA", colour_by = "All|MHb.1.2"
) + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(NA, 4), oob = scales::squish) +
  ggtitle('Human MHb.1.2 markers in Hashikawa')

plotReducedDim(
  ortho_present_hashikawa_sce_sub,
  dimred = "TSNE.PCA", colour_by = "All|MHb.2"
) + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(NA, 4), oob = scales::squish) +
  ggtitle('Human MHb.2 markers in Hashikawa')

plotReducedDim(
  ortho_present_hashikawa_sce_sub,
  dimred = "TSNE.PCA", colour_by = "All|LHb.1.3.4"
) + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(NA, 4), oob = scales::squish) +
  ggtitle('Human LHb.1.3.4 markers in Hashikawa')

plotReducedDim(
  ortho_present_hashikawa_sce_sub,
  dimred = "TSNE.PCA", colour_by = "All|LHb.2.7"
) + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(NA, 4), oob = scales::squish) +
  ggtitle('Human LHb.2.7 markers in Hashikawa')

plotReducedDim(
  ortho_present_hashikawa_sce_sub,
  dimred = "TSNE.PCA", colour_by = "All|LHb.4"
) + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(NA, 4), oob = scales::squish) +
  ggtitle('Human LHb.4 markers in Hashikawa')






