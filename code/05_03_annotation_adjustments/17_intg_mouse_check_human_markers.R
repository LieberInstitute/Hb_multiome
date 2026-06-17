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


both_mouse_sce = cbind(all_mouse_sce, hashikawa_sce_sub)

#Filter out the outlier and Neuron cells
both_mouse_sce = both_mouse_sce[, !both_mouse_sce$grouped_annot %in% c('outliers', 'Neuron')]


both_mouse_sce = logNormCounts(both_mouse_sce,name = "logcounts")

#And standard HVG
dec.both_mouse_sce <- modelGeneVar(both_mouse_sce, assay.type = "logcounts")
hvgs <- getTopHVGs(dec.both_mouse_sce, n = 2000)


both_mouse_sce <- runPCA(both_mouse_sce,
    subset_row = hvgs,
    ncomponents = 30,
    name = "PCA"
)




## Run Harmony
message("Running Harmony - ", Sys.time())
both_mouse_sce <- RunHarmony(both_mouse_sce, group.by.vars = c('batch_variable_2'), verbose = TRUE)


#Visualize before and after correction
message("Running TSNE - ", Sys.time())
both_mouse_sce <- runTSNE(both_mouse_sce, dimred = "HARMONY", name = "TSNE.HARMONY")
colnames(reducedDim(both_mouse_sce, "TSNE.HARMONY")) <- c("TSNE1", "TSNE2")

message("Running UMAP - ", Sys.time())
both_mouse_sce <- runUMAP(both_mouse_sce, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(both_mouse_sce, "UMAP.HARMONY")) <- c("UMAP1", "UMAP2")


message("Running TSNE - ", Sys.time())
both_mouse_sce <- runTSNE(both_mouse_sce, dimred = "PCA", name = "TSNE.PCA")
colnames(reducedDim(both_mouse_sce, "TSNE.PCA")) <- c("TSNE1", "TSNE2")

message("Running UMAP - ", Sys.time())
both_mouse_sce <- runUMAP(both_mouse_sce, dimred = "PCA", name = "UMAP.PCA")
colnames(reducedDim(both_mouse_sce, "UMAP.PCA")) <- c("UMAP1", "UMAP2")


both_mouse_sce$grouped_annot[both_mouse_sce$grouped_annot %in% c()]


plotReducedDim(
  both_mouse_sce,
  dimred = "PCA", colour_by = "batch_variable_2"
)

plotReducedDim(
  both_mouse_sce,
  dimred = "HARMONY", colour_by = "batch_variable"
)


plotReducedDim(
  both_mouse_sce,
  dimred = "TSNE.HARMONY", colour_by = "batch_variable"
)

plotReducedDim(
  both_mouse_sce,
  dimred = "TSNE.HARMONY", colour_by = "grouped_annot"
)


#And now visualize the top human markers in this integrated data.







