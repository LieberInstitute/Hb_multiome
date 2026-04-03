#Comparing the habenula metaclusters across zebrafish and mouse with MetaNeighbor
#Might be difficult due to the zebrafish genome duplication, a lot of informative genes have one-many matches between mouse and zebrafish
#Those one-many genes will get filtered out with the traditional approach of just keeping the one-to-ones

library(SingleCellExperiment)
library(Seurat)
library(MetaNeighbor)
library(dplyr)
library(ggplot2)
library(here)

here::here()


#Paths to Wallace annotated data
wallace_02_data_path = here('processed-data','98_external_Hb_comparisons','02_qc_and_clust_wallace_2019')
wallace_03_data_path = here('processed-data','98_external_Hb_comparisons','03_metaMarkers_wallace_2019')

#Paths to Hashikawa annotated data
hashikawa_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/09_cross_species_analysis/Hashikawa_data'

#Path to the Hashikawa metadata I want
hashikawa_data_path = here('processed-data','98_external_Hb_comparisons','04_wallace_hashikawa_mouse')

#Paths to Pendey annotated data
pendey_05_data_path = here('processed-data','98_external_Hb_comparisons','05_initial_qc_pandey_2018')
pendey_06_data_path = here('processed-data', '98_external_Hb_comparisons', '06_metaMarkers_pandey_2018')


#Path to orthologs
path_to_orthologs = here('processed-data','98_external_Hb_comparisons','human_mouse_zebrafish_orthologs.txt.gz')


new_data_path = here('processed-data','98_external_Hb_comparisons','07_cross_species_hab')
plot_path = here('plots','98_external_Hb_comparisons','07_cross_species_hab')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)



#Load up the full SCE objects, contains the metacluster annotations
all_mouse_sce = readRDS(paste0(wallace_02_data_path, '/all_donor_sce_with_denovo_clusters.rds'))

all_zebF_sce = readRDS(paste0(pendey_05_data_path, '/all_donor_sce_with_denovo_clusters.rds'))

#Hashikawa data
#Loads as sce_mouse_sub, list of sce objects
load(paste0(hashikawa_path, '/sce_mouse_habenula.Rdata'))
hashikawa_sce_sub = sce_mouse_sub$neuron #Using just the neuron subset
rm(sce_mouse_sub)

full_hashikawa_metadata = readRDS(paste0(hashikawa_data_path, '/hashikawa_mouse_neuron_metaclust_celltype_annot_metadata.rds'))
colnames(full_hashikawa_metadata )
colData(hashikawa_sce_sub) = full_hashikawa_metadata 
#Gene symbols as the rownames
rownames(hashikawa_sce_sub) = rowData(hashikawa_sce_sub)$Symbol
#Add CPM
assay(hashikawa_sce_sub, "cpm") = MetaMarkers::convert_to_cpm(assay(hashikawa_sce_sub, "counts"))



#Filter for the present 1-to-1 orthologs, might as well use the set of orthologs between human, mouse, and zebrafish
hu_mu_zf_ortholog_df = data.table::fread(path_to_orthologs)
table(hu_mu_zf_ortholog_df$`Mouse homology type`)

#Get the 1to1 orthologs, 9435 in total
hu_mu_zf_ortholog_df <- hu_mu_zf_ortholog_df %>%
  filter(
    `Mouse homology type` == 'ortholog_one2one' &
      `Zebrafish homology type` == 'ortholog_one2one'
  ) %>%
  filter(!duplicated(`Gene name`))
dim(hu_mu_zf_ortholog_df)

#Swap in the human gene names for each specie
index = match(
  rownames(all_zebF_sce),
  toupper(hu_mu_zf_ortholog_df$`Zebrafish gene name`)
)
rownames(all_zebF_sce) = hu_mu_zf_ortholog_df$`Gene name`[index]

index = match(rownames(all_mouse_sce), hu_mu_zf_ortholog_df$`Mouse gene name`)
rownames(all_mouse_sce) = hu_mu_zf_ortholog_df$`Gene name`[index]

index = match(rownames(hashikawa_sce_sub), hu_mu_zf_ortholog_df$`Mouse gene name`)
rownames(hashikawa_sce_sub) = hu_mu_zf_ortholog_df$`Gene name`[index]


#Get the set of genes that are present in both datasets, and are 1-to-1 orthologs
all_present_genes = intersect(intersect(rownames(all_mouse_sce), rownames(all_zebF_sce)), rownames(hashikawa_sce_sub))
all_present_genes = all_present_genes[!is.na(all_present_genes)]
length(all_present_genes)

#Filter the SCE objects for the present genes
#7415 genes total
all_zebF_sce = all_zebF_sce[rownames(all_zebF_sce) %in% all_present_genes, ]
all_mouse_sce = all_mouse_sce[rownames(all_mouse_sce) %in% all_present_genes, ]
hashikawa_sce_sub = hashikawa_sce_sub[rownames(hashikawa_sce_sub) %in% all_present_genes, ]

dim(all_zebF_sce)
dim(all_mouse_sce)
dim(hashikawa_sce_sub)
#Sanity check
table(rownames(all_zebF_sce) %in% rownames(all_mouse_sce))
table(rownames(all_zebF_sce) %in% rownames(hashikawa_sce_sub))


#add in the annotated metdata
current_mouse_metadata = readRDS(paste0(
  wallace_03_data_path,
  '/wallace_mouse_metaclust_celltype_annot_metadata.rds'
))
colData(all_mouse_sce) = S4Vectors::DataFrame(current_mouse_metadata)

#Now do the same for the zebrafish data
current_zebF_metadata = readRDS(paste0(
  pendey_06_data_path,
  '/pandey_zebrafish_metaclust_celltype_annot_metadata.rds'
))
colData(all_zebF_sce) = S4Vectors::DataFrame(current_zebF_metadata)


#Just focus on the neuronal populations
all_mouse_sce = all_mouse_sce[, !all_mouse_sce$meta_clust_celltype_annot %in% c('Astrocytes','Polydendrocytes','Differentiating Oligodendrocytes','Oligodendrocytes',
'Pericytes','Fibroblasts','Endothelial','Macrophages', 'Microglia')]
table(all_mouse_sce$meta_clust_celltype_annot)


all_zebF_sce = all_zebF_sce[, !all_zebF_sce$meta_clust_celltype_annot %in% c('non_neuronal')]
table(all_zebF_sce$meta_clust_celltype_annot)


#Split and then merge
table(all_mouse_sce$study_id)

table(all_zebF_sce$study_id)

studies <- unique(all_mouse_sce$study_id)
mouse_sce_list <- lapply(studies, function(study) {
  all_mouse_sce[, all_mouse_sce$study_id == study]
})
names(mouse_sce_list) <- studies


studies <- unique(all_zebF_sce$study_id)
zebF_sce_list <- lapply(studies, function(study) {
  all_zebF_sce[, all_zebF_sce$study_id == study]
})
names(zebF_sce_list) <- studies

#Split Hashikawa by stim or control
studies <- unique(hashikawa_sce_sub$stim)
hashikawa_sce_sub_list <- lapply(studies, function(study) {
  hashikawa_sce_sub[, hashikawa_sce_sub$stim == study]
})
names(hashikawa_sce_sub_list) <- studies


########################
#MetaNeighbor
#########################

#Get single SCE object
all_donor_sce = mergeSCE(c(mouse_sce_list, zebF_sce_list, hashikawa_sce_sub_list))
View(as.data.frame(colData(all_donor_sce)))

#Ignore the outlier cells
all_donor_sce = all_donor_sce[,
  all_donor_sce$meta_clust_celltype_annot != 'outliers'
]



#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(
  dat = all_donor_sce,
  min_recurrence = 2,
  exp_labels = all_donor_sce$study_id
)
length(global_hvgs)
keep_global_hvgs = global_hvgs[1:2000]


MN_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$meta_clust_celltype_annot,
  fast_version = TRUE
)

#Plot allby-all AUROC heatmap
plotHeatmap(
  MN_aurocs,
  cex = .5
)
title("MetaNeighbor Mouse and Zebrafish Habenula: 2000 HVGs non-summed paralogs")

pdf(paste0(plot_path, '/AllvsAll_MN_mouse_zeb_neurons.pdf'), width = 10, height = 8)
plotHeatmap(
  MN_aurocs,
  cex = .5
)
title("MetaNeighbor Mouse and Zebrafish Habenula: 2000 HVGs non-summed paralogs")
dev.off()




#And the best versus next approach
MN_best_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$meta_clust_celltype_annot,
  fast_version = TRUE,
  one_vs_best = TRUE,
  symmetric_output = FALSE
)

#Plot best_vs_next AUROC heatmap
plotHeatmap(
  MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Mouse and Zebrafish Habenula: 2000 HVGs non-summed paralogs")

pdf(paste0(plot_path, '/BvsNext_MN_mouse_zeb_neurons.pdf'), width = 10, height = 8)
plotHeatmap(
  MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Mouse and Zebrafish Habenula: 2000 HVGs non-summed paralogs")
dev.off()


cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$study_id, all_donor_sce$meta_clust_celltype_annot, size_factor = 3)

pdf(paste0(plot_path, '/MN_cluster_graph.pdf'), width = 10, height = 8)
cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .7)
plotClusterGraph(cluster_graph, all_donor_sce$study_id, all_donor_sce$meta_clust_celltype_annot, size_factor = 3)
dev.off()


#Get the metaclusters from the best vs next results, add those annotations to the full SCE object
mclusters = extractMetaClusters(MN_best_aurocs, threshold = .9)
mclusters



##################################
#
#Try this again but with a little work on the zebrafish data
#A lot of the informative marker genes used were one-many orthologs between mouse and zebrafish
#Due to the teleost genome duplication
#Since there's no easy way to tell which A or B copy is the actual best comparison to the mouse gene
#Let's try just summing the counts between A and B paralogs in zebrafish, then matching to the 1-to-1 ortholog between human and mouse

#Reload data for the full gene set
all_zebF_sce = readRDS(paste0(
  pendey_05_data_path,
  '/all_donor_sce_with_denovo_clusters.rds'
))


#Get the human gene name, using the one-to-many orthologs in zebrafish
hu_mu_zf_ortholog_df = data.table::fread(path_to_orthologs)
hu_mu_zf_ortholog_df <- hu_mu_zf_ortholog_df %>%
  filter(
    `Mouse homology type` == 'ortholog_one2one' &
      `Zebrafish homology type` %in% c('ortholog_one2one', 'ortholog_one2many')
  ) %>%
  filter(!duplicated(`Zebrafish gene name`)) #Want to keep the human duplicate gene names
dim(hu_mu_zf_ortholog_df)
View(hu_mu_zf_ortholog_df)

#Turn the zebrafish counts matrix into a dataframe so I can use simple dplyr functions

zeb_counts_df = as.data.frame(assay(all_zebF_sce))
zeb_counts_df$zeb_gene_name = rownames(zeb_counts_df)
#Add the human gene name 
index = match(
  zeb_counts_df$zeb_gene_name,
  toupper(hu_mu_zf_ortholog_df$`Zebrafish gene name`)
)
zeb_counts_df$human_gene_name = hu_mu_zf_ortholog_df$`Gene name`[index]

#Check that the grouping is working
#88 human genes have more than 2 sebrafish genes mapping to them, minor number but summing all those might still be weird
zeb_counts_df %>% filter(!is.na(human_gene_name)) %>% 
  group_by(human_gene_name) %>% summarise(count = n()) %>% filter(count > 2)

#Get a test matrix to double check the grouping and summing is working
test_counts = zeb_counts_df %>%
  filter(zeb_gene_name %in% c('SLC17A7A','SLC17A7B')) %>% 
  select(-zeb_gene_name, -human_gene_name)
test_counts = test_counts[ ,colSums(test_counts) > 1][1:20]  #20 cells with non-zero counts
example_colsums_test = colSums(test_counts)
test_counts$zeb_gene_name = c('SLC17A7A','SLC17A7B')
test_counts$human_gene_name = c('SLC17A7','SLC17A7')


group_sum_example = test_counts %>%
  filter(!is.na(human_gene_name)) %>%
  select(-zeb_gene_name) %>%
  group_by(human_gene_name) %>%
  summarise(across(everything(), sum)) %>%
  tibble::column_to_rownames("human_gene_name") %>%
  as.matrix()


table(group_sum_example == example_colsums_test) #All the same, seems to be working


#Seems to work, will probably take a little while to run on the full dataset
zeb_sum_counts_mat = zeb_counts_df %>%
  filter(!is.na(human_gene_name)) %>%
  select(-zeb_gene_name) %>%
  group_by(human_gene_name) %>%
  summarise(across(everything(), sum)) %>%
  tibble::column_to_rownames("human_gene_name") %>%
  as.matrix()

#Drop the first row, empty string gene name
zeb_sum_counts_mat = zeb_sum_counts_mat[2:nrow(zeb_sum_counts_mat), ]


dim(zeb_counts_df )
dim(zeb_sum_counts_mat)


#Make a new seurat experiment with the summed counts, and the annotated metadata
current_zebF_metadata = readRDS(paste0(
  pendey_06_data_path,
  '/pandey_zebrafish_metaclust_celltype_annot_metadata.rds'
))



zeb_seurat <- CreateSeuratObject(counts = zeb_sum_counts_mat, meta.data = current_zebF_metadata)
View(zeb_seurat[[]])

#Full processing just to see the umap, if wildly different
 zeb_seurat <- NormalizeData(zeb_seurat , normalization.method = "RC", scale.factor = 1e6)
#HVGs
zeb_seurat <- FindVariableFeatures(zeb_seurat , selection.method = "vst", nfeatures = 2000)

#Scale data
all.genes <- rownames(zeb_seurat )
zeb_seurat <- ScaleData(zeb_seurat , features = all.genes)

#PCA
zeb_seurat  <- RunPCA(zeb_seurat , features = VariableFeatures(object = zeb_seurat ))
DimPlot(zeb_seurat , reduction = "pca") + NoLegend()

#UMAP
zeb_seurat  <- RunUMAP(zeb_seurat  , dims = 1:20)

#Doesnt look insane
p1 = DimPlot(zeb_seurat , reduction = "umap", group.by = 'meta_clust_celltype_annot', label = TRUE) + 
  ggtitle('Pendey 2018: MetaCluster annotations, summed paralogs')
p2 = DimPlot(zeb_seurat , reduction = "umap", group.by = 'study_id', label = TRUE) + 
  ggtitle('Pendey 2018: MetaCluster annotations, summed paralogs')

p1
p2

ggsave(p1, filename = 'pandey_zebrafish_hab_summed_paralog_annot_umap.pdf', path = plot_path,
device = 'pdf', width = 8, height = 7)
ggsave(p2, filename = 'pandey_zebrafish_hab_summed_paralog_batch_umap.pdf', path = plot_path,
device = 'pdf', width = 8, height = 7)


#Convert to SingleCellExperiment for MetaNeighbor
#And also filter the mouse data to match, the new zebrafish dataset is already in human gene names
zeb_sce = as.SingleCellExperiment(zeb_seurat)
names(assays(zeb_sce)) = c('counts','cpm','scaledata')
zeb_sce

#Save the paralog summed version of the zebrafish data
saveRDS(zeb_sce, file = paste0(new_data_path, '/adult_zebrafish_summed_paralogs.rds'))




#Load up the full SCE objects, contains the metacluster annotations
all_mouse_sce = readRDS(paste0(
  wallace_02_data_path,
  '/all_donor_sce_with_denovo_clusters.rds'
))

##Hashikawa data
load(paste0(hashikawa_path, '/sce_mouse_habenula.Rdata'))
hashikawa_sce_sub = sce_mouse_sub$neuron #Using just the neuron subset
rm(sce_mouse_sub)

full_hashikawa_metadata = readRDS(paste0(hashikawa_data_path, '/hashikawa_mouse_neuron_metaclust_celltype_annot_metadata.rds'))
colnames(full_hashikawa_metadata )
colData(hashikawa_sce_sub) = full_hashikawa_metadata 
#Gene symbols as the rownames
rownames(hashikawa_sce_sub) = rowData(hashikawa_sce_sub)$Symbol
#Add CPM
assay(hashikawa_sce_sub, "cpm") = MetaMarkers::convert_to_cpm(assay(hashikawa_sce_sub, "counts"))



#For the mouse to human genes 
hu_mu_zf_ortholog_df = data.table::fread(path_to_orthologs)
table(hu_mu_zf_ortholog_df$`Mouse homology type`)

#Get the 1to1 orthologs
hu_mu_zf_ortholog_df <- hu_mu_zf_ortholog_df %>%
  filter(
    `Mouse homology type` == 'ortholog_one2one'
  ) %>%
  filter(!duplicated(`Gene name`))
dim(hu_mu_zf_ortholog_df)

index = match(rownames(all_mouse_sce), hu_mu_zf_ortholog_df$`Mouse gene name`)
rownames(all_mouse_sce) = hu_mu_zf_ortholog_df$`Gene name`[index]

index = match(rownames(hashikawa_sce_sub), hu_mu_zf_ortholog_df$`Mouse gene name`)
rownames(hashikawa_sce_sub) = hu_mu_zf_ortholog_df$`Gene name`[index]



#Get the set of genes that are present in both datasets, and are 1-to-1 orthologs
all_present_genes = intersect(intersect(rownames(all_mouse_sce), rownames(zeb_sce)), rownames(hashikawa_sce_sub))
all_present_genes = all_present_genes[!is.na(all_present_genes)]
length(all_present_genes)
#Filter the SCE objects for the present genes
#9532 genes total
zeb_sce = zeb_sce[rownames(zeb_sce) %in% all_present_genes, ]
all_mouse_sce = all_mouse_sce[rownames(all_mouse_sce) %in% all_present_genes, ]
hashikawa_sce_sub = hashikawa_sce_sub[rownames(hashikawa_sce_sub) %in% all_present_genes, ]


dim(zeb_sce)
dim(all_mouse_sce)
dim(hashikawa_sce_sub)
#Sanity check
table(rownames(zeb_sce) %in% rownames(all_mouse_sce))
table(rownames(zeb_sce) %in% rownames(hashikawa_sce_sub))


#Add the mouse metadata
current_mouse_metadata = readRDS(paste0(
  wallace_03_data_path,
  '/wallace_mouse_metaclust_celltype_annot_metadata.rds'
))
colData(all_mouse_sce) = S4Vectors::DataFrame(current_mouse_metadata)


#Just focus on the neuronal populations
all_mouse_sce = all_mouse_sce[, !all_mouse_sce$meta_clust_celltype_annot %in% c('Astrocytes','Polydendrocytes','Differentiating Oligodendrocytes','Oligodendrocytes',
'Pericytes','Fibroblasts','Endothelial','Macrophages', 'Microglia')]
table(all_mouse_sce$meta_clust_celltype_annot)


zeb_sce = zeb_sce[, !zeb_sce$meta_clust_celltype_annot %in% c('non_neuronal')]
table(zeb_sce$meta_clust_celltype_annot)



studies <- unique(all_mouse_sce$study_id)
mouse_sce_list <- lapply(studies, function(study) {
  all_mouse_sce[, all_mouse_sce$study_id == study]
})
names(mouse_sce_list) <- studies


studies <- unique(zeb_sce$study_id)
zebF_sce_list <- lapply(studies, function(study) {
  zeb_sce[, zeb_sce$study_id == study]
})
names(zebF_sce_list) <- studies


#Split Hashikawa by stim or control
studies <- unique(hashikawa_sce_sub$stim)
hashikawa_sce_sub_list <- lapply(studies, function(study) {
  hashikawa_sce_sub[, hashikawa_sce_sub$stim == study]
})
names(hashikawa_sce_sub_list) <- studies

########################
#MetaNeighbor
#########################

#Get single SCE object
all_donor_sce = mergeSCE(c(mouse_sce_list, zebF_sce_list, hashikawa_sce_sub_list))
View(as.data.frame(colData(all_donor_sce)))

#Ignore the outlier cells
all_donor_sce = all_donor_sce[,
  all_donor_sce$meta_clust_celltype_annot != 'outliers'
]

#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(
  dat = all_donor_sce,
  min_recurrence = 2,
  exp_labels = all_donor_sce$study_id
)
length(global_hvgs)
keep_global_hvgs = global_hvgs[1:2000]

keep_global_hvgs[grepl('GAD', keep_global_hvgs)]
keep_global_hvgs[grepl('SLC17A', keep_global_hvgs)]


paralogSummed_MN_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$meta_clust_celltype_annot,
  fast_version = TRUE
)

#Plot allby-all AUROC heatmap
plotHeatmap(
  paralogSummed_MN_aurocs,
  cex = .5
)
title("MetaNeighbor Mouse and Zebrafish Habenula: 2000 HVGs summed paralogs")

pdf(paste0(plot_path, '/paralogSummed_AllvsAll_MN_mouse_zeb_neurons.pdf'), width = 10, height = 8)
plotHeatmap(
  paralogSummed_MN_aurocs,
  cex = .5
)
title("MetaNeighbor Mouse and Zebrafish Habenula: 2000 HVGs summed paralogs")
dev.off()


#And the best versus next approach
paralogSummed_MN_best_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$meta_clust_celltype_annot,
  fast_version = TRUE,
  one_vs_best = TRUE,
  symmetric_output = FALSE
)

#Plot best_vs_next AUROC heatmap
plotHeatmap(
  paralogSummed_MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Mouse and Zebrafish Habenula: 2000 HVGs summed paralogs")


pdf(paste0(plot_path, '/paralogSummed_BvsNext_MN_mouse_zeb_neurons.pdf'), width = 10, height = 8)
plotHeatmap(
  paralogSummed_MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Mouse and Zebrafish Habenula: 2000 HVGs summed paralogs")
dev.off()



paralog_cluster_graph = makeClusterGraph(paralogSummed_MN_best_aurocs, low_threshold = .7)
plotClusterGraph(paralog_cluster_graph, all_donor_sce$study_id, all_donor_sce$meta_clust_celltype_annot, size_factor = 3)

pdf(paste0(plot_path, '/paralogSummed_MN_cluster_graph.pdf'), width = 10, height = 8)
plotClusterGraph(paralog_cluster_graph, all_donor_sce$study_id, all_donor_sce$meta_clust_celltype_annot, size_factor = 3)
dev.off()




#Get the metaclusters from the best vs next results, add those annotations to the full SCE object
mclusters = extractMetaClusters(paralogSummed_MN_best_aurocs, threshold = .5)
mclusters



