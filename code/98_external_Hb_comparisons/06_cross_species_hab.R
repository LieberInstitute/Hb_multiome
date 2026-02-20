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
wallace_02_data_path = here('processed-data', '98_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')
wallace_03_data_path = here('processed-data', '98_external_Hb_comparisons', '03_metaMarkers_wallace_2019')

#Paths to Pendey annotated data
pendey_04_data_path = here('processed-data', '98_external_Hb_comparisons', '04_initial_qc_pandey_2018')
pendey_05_data_path = here('processed-data', '98_external_Hb_comparisons', '05_metaMarkers_pandey_2018')

#Path to orthologs
path_to_orthologs = here('processed-data', '98_external_Hb_comparisons', 'human_mouse_zebrafish_orthologs.txt.gz')



#Load up the full SCE objects, contains the metacluster annotations
all_mouse_sce = readRDS(paste0(wallace_02_data_path, '/all_donor_sce_with_denovo_clusters.rds'))
all_zebF_sce = readRDS(paste0(pendey_04_data_path, '/all_donor_sce_with_denovo_clusters.rds'))


#Filter for the present 1-to-1 orthologs, might as well use the set of orthologs between human, mouse, and zebrafish
hu_mu_zf_ortholog_df = data.table::fread(path_to_orthologs)
table(hu_mu_zf_ortholog_df$`Mouse homology type`)

#Get the 1to1 orthologs, 9435 in total
hu_mu_zf_ortholog_df <- hu_mu_zf_ortholog_df %>% 
  filter(`Mouse homology type` == 'ortholog_one2one' & `Zebrafish homology type` == 'ortholog_one2one') %>%
  filter(!duplicated(`Gene name`))
dim(hu_mu_zf_ortholog_df)

#Swap in the human gene names for each specie
index = match(rownames(all_zebF_sce), toupper(hu_mu_zf_ortholog_df$`Zebrafish gene name`))
rownames(all_zebF_sce) = hu_mu_zf_ortholog_df$`Gene name`[index]

index = match(rownames(all_mouse_sce), hu_mu_zf_ortholog_df$`Mouse gene name`)
rownames(all_mouse_sce) = hu_mu_zf_ortholog_df$`Gene name`[index]

#Get the set of genes that are present in both datasets, and are 1-to-1 orthologs
all_present_genes = intersect(rownames(all_mouse_sce), rownames(all_zebF_sce))
all_present_genes = all_present_genes[!is.na(all_present_genes)]

#Filter the SCE objects for the present genes
#8142 genes total
all_zebF_sce = all_zebF_sce[rownames(all_zebF_sce) %in% all_present_genes, ]
all_mouse_sce = all_mouse_sce[rownames(all_mouse_sce) %in% all_present_genes, ]

dim(all_zebF_sce)
dim(all_mouse_sce)
#Sanity check
table(rownames(all_zebF_sce) %in% rownames(all_mouse_sce))



#Convert to seurat and then add in the annotated metdata
current_mouse_metadata = readRDS(paste0(wallace_03_data_path, '/wallace_mouse_metaclust_celltype_annot_metadata.rds'))
colData(all_mouse_sce) = S4Vectors::DataFrame(current_mouse_metadata)

#Now do the same for the zebrafish data
current_zebF_metadata = readRDS(paste0(pendey_05_data_path, '/pandey_zebrafish_metaclust_celltype_annot_metadata.rds'))
colData(all_zebF_sce) = S4Vectors::DataFrame(current_zebF_metadata)



#Seems silly to split than merge for MetaNeighbor, but I want to make sure the cluster labels are in the same column for both datasets
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



########################
#MetaNeighbor
#########################

#Get single SCE object
all_donor_sce = mergeSCE(c(mouse_sce_list, zebF_sce_list))
View(as.data.frame(colData(all_donor_sce)))

#Ignore the outlier cells
all_donor_sce = all_donor_sce[, all_donor_sce$meta_clust_celltype_annot != 'outliers']
#Ignore larva?
all_donor_sce = all_donor_sce[, all_donor_sce$study_id != 'Larva']
table(all_donor_sce$study_id)



#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(dat = all_donor_sce, min_recurrence = 2, exp_labels = all_donor_sce$study_id)
length(global_hvgs)
keep_global_hvgs = global_hvgs[1:3000]


MN_aurocs = MetaNeighborUS(var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$meta_clust_celltype_annot,
  fast_version = TRUE)

#Plot allby-all AUROC heatmap
plotHeatmap(MN_aurocs, 
  show_dendro = TRUE, 
  show_labels = TRUE, 
  cex = .5,
  title = "MetaNeighbor AUROCs for Mouse and Zebrafish Habenula")


#And the best versus next approach

MN_best_aurocs = MetaNeighborUS(var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$study_id,
  cell_type = all_donor_sce$meta_clust_celltype_annot,
  fast_version = TRUE,
one_vs_best = TRUE, symmetric_output = FALSE)

#Plot best_vs_next AUROC heatmap
plotHeatmap(MN_best_aurocs, 
  show_dendro = TRUE, 
  show_labels = TRUE, 
  cex = .5,
  title = "MetaNeighbor best_vs_next AUROCs Mouse and Zebrafish Habenula")

#Get the metaclusters from the best vs next results, add those annotations to the full SCE object
mclusters = extractMetaClusters(MN_best_aurocs, threshold = .9)
mclusters






