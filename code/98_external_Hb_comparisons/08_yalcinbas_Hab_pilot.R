#Loading up and seeing what is currently present in the published Yalcinbas Habenula single cell data
#The single cell data is all neurotypical donors


library(SingleCellExperiment)
library(Seurat)
library(MetaNeighbor)
library(dplyr)
library(ggplot2)
library(here)

here::here()


#Path to the Yalcinbas pilot data
#Going with the official_final_sce.RDATA
yalcinbas_path = '/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/sce_objects'
list.files(yalcinbas_path)

#Nonhuman data paths
zeb_data_path = here('processed-data', '98_external_Hb_comparisons','07_cross_species_hab')
mouse_data_path = here('processed-data', '98_external_Hb_comparisons', '02_qc_and_clust_wallace_2019')
wallace_03_data_path = here('processed-data', '98_external_Hb_comparisons', '03_metaMarkers_wallace_2019')

#Path to orthologs
path_to_orthologs = here('processed-data', '98_external_Hb_comparisons', 'human_mouse_zebrafish_orthologs.txt.gz')


#Path for new data generated
new_data_path = here('processed-data', '98_external_Hb_comparisons','08_yalcinbas_Hab_pilot')
#Path to plot directory
plot_path = here('plots', '98_external_Hb_comparisons', '08_yalcinbas_Hab_pilot')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)

#Loads as an object labeled 'sce'
#16437 cells
load(paste0(yalcinbas_path, '/official_final_sce.RDATA'))
yalcinbas_sce = sce
rm(sce)

colnames(colData(yalcinbas_sce))
table(yalcinbas_sce$final_Annotations)
table(yalcinbas_sce$broad_Annotations)

table(yalcinbas_sce$Sample, yalcinbas_sce$final_Annotations )

yalcinbas_sce


#Just try the whole dataset as a sample, can break it up by donor later, but there are only 4 donors with decent annotated habenula cells

#Load up the zebrafish and mouse data, filter down to the present genes
#Zebrafish data, using the summed paralog version
zeb_sce = readRDS(file = paste0(zeb_data_path, '/adult_zebrafish_summed_paralogs.rds'))

#Wallace mouse data
all_mouse_sce = readRDS(paste0(mouse_data_path, '/all_donor_sce_with_denovo_clusters.rds'))



hu_mu_zf_ortholog_df = data.table::fread(path_to_orthologs)
table(hu_mu_zf_ortholog_df$`Mouse homology type`)

#Get the 1to1 orthologs for the mouse
hu_mu_zf_ortholog_df <- hu_mu_zf_ortholog_df %>%
  filter(
    `Mouse homology type` == 'ortholog_one2one'
  ) %>%
  filter(!duplicated(`Gene name`))
dim(hu_mu_zf_ortholog_df)

index = match(rownames(all_mouse_sce), hu_mu_zf_ortholog_df$`Mouse gene name`)
rownames(all_mouse_sce) = hu_mu_zf_ortholog_df$`Gene name`[index]

#Get the set of genes that are present in all datasets, and are 1-to-1 orthologs
#10431 across all three, which is pretty good
all_present_genes = Reduce(intersect, list(rownames(all_mouse_sce), rownames(zeb_sce), rownames(yalcinbas_sce)))
all_present_genes = all_present_genes[!is.na(all_present_genes)]
length(all_present_genes)


#Filter the SCE objects for the present genes
#10431 genes total
zeb_sce = zeb_sce[rownames(zeb_sce) %in% all_present_genes, ]
all_mouse_sce = all_mouse_sce[rownames(all_mouse_sce) %in% all_present_genes, ]
yalcinbas_sce = yalcinbas_sce[rownames(yalcinbas_sce) %in% all_present_genes, ]

dim(zeb_sce)
dim(all_mouse_sce)
dim(yalcinbas_sce)
#Sanity check
table(rownames(zeb_sce) %in% rownames(all_mouse_sce))
table(rownames(zeb_sce) %in% rownames(yalcinbas_sce))


#Add the mouse metadata
current_mouse_metadata = readRDS(paste0(wallace_03_data_path, '/wallace_mouse_metaclust_celltype_annot_metadata.rds'))
colData(all_mouse_sce) = S4Vectors::DataFrame(current_mouse_metadata)


#Rename and add metadata so the columns match across the datasets
#Rename the meta_clust_celltype_annot columns in the zebrafish and mouse data to final_Annotations to match Yalcinbas
zeb_sce$final_Annotations = zeb_sce$meta_clust_celltype_annot
all_mouse_sce$final_Annotations = all_mouse_sce$meta_clust_celltype_annot

zeb_sce$species = 'Zebrafish'
all_mouse_sce$species = 'Mouse'
yalcinbas_sce$species = 'Human'


#Split the mouse and zebrafish SCE objects into lists of SCE objects for each donor, for MetaNeighbor
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


#Get single SCE object
all_donor_sce = mergeSCE(c(mouse_sce_list, zebF_sce_list, list(Human = yalcinbas_sce)))
table(all_donor_sce$study_id)
table(all_donor_sce$species)
table(all_donor_sce$final_Annotations)

#Ignore the outlier cells
all_donor_sce = all_donor_sce[, all_donor_sce$final_Annotations!= 'outliers']


#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(
  dat = all_donor_sce,
  min_recurrence = 2,
  exp_labels = all_donor_sce$species
)
length(global_hvgs)
#keep_global_hvgs = global_hvgs[1:2000]
keep_global_hvgs = global_hvgs

MN_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$species,
  cell_type = all_donor_sce$final_Annotations,
  fast_version = TRUE
)

#Plot allby-all AUROC heatmap
plotHeatmap(
  MN_aurocs,
  cex = .5
)
title("MetaNeighbor Human, Mouse, ZebF Habenula: 2000 HVGs")

#And the best versus next approach

MN_best_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = all_donor_sce,
  study_id = all_donor_sce$species,
  cell_type = all_donor_sce$final_Annotations,
  fast_version = TRUE,
  one_vs_best = TRUE,
  symmetric_output = FALSE
)

#Plot best_vs_next AUROC heatmap
plotHeatmap(
  MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Human, Mouse, ZebF Habenula: 2000 HVGs")

cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .3)
plotClusterGraph(cluster_graph, all_donor_sce$species, all_donor_sce$final_Annotations, size_factor = 3)





#And what about just mouse and human
#Reload and get the overlap of orthologs just for mouse and human

#Get single SCE object
human_mouse_sce = mergeSCE(c(mouse_sce_list, list(Human = yalcinbas_sce)))
table(human_mouse_sce$study_id)
table(human_mouse_sce$species)
table(human_mouse_sce$final_Annotations)

#Ignore the outlier cells
human_mouse_sce = human_mouse_sce[, human_mouse_sce$final_Annotations!= 'outliers']


#Get highly variable genes, this time highly variable genes across the donor datasets, sticking with 2000
global_hvgs = variableGenes(
  dat = human_mouse_sce,
  min_recurrence = 2,
  exp_labels = human_mouse_sce$species
)
length(global_hvgs)
keep_global_hvgs = global_hvgs[1:800]


MN_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = human_mouse_sce,
  study_id = human_mouse_sce$species,
  cell_type = human_mouse_sce$final_Annotations,
  fast_version = TRUE
)

#Plot allby-all AUROC heatmap
plotHeatmap(
  MN_aurocs,
  cex = .5
)
title("MetaNeighbor Human, Mouse Habenula: 2000 HVGs")

#And the best versus next approach

MN_best_aurocs = MetaNeighborUS(
  var_genes = keep_global_hvgs,
  dat = human_mouse_sce,
  study_id = human_mouse_sce$species,
  cell_type = human_mouse_sce$final_Annotations,
  fast_version = TRUE,
  one_vs_best = TRUE,
  symmetric_output = FALSE
)

#Plot best_vs_next AUROC heatmap
plotHeatmap(
  MN_best_aurocs,
  cex = .5
)
title("MetaNeighbor BvsNext Human, Mouse Habenula: 2000 HVGs")



#From the most obvious matches between mouse and human, I would predict that 
#Human MHb.2 would be cholinergic
#Human MHb.3 is likely cholinergic

#Human MHb.1 would be substnace P

#Human LHb.6 is a closer match to human medial clusters than lateral ones

#Mouse LHb_1 has a weak match in human LHb.3 and .4, 
#Mouse LHb_2 has a weak match in human LHb.2

#Get the metaclusters from the best vs next results, add those annotations to the full SCE object
mclusters = extractMetaClusters(MN_best_aurocs, threshold = .5)
mclusters

cluster_graph = makeClusterGraph(MN_best_aurocs, low_threshold = .3)
plotClusterGraph(cluster_graph, human_mouse_sce$species, human_mouse_sce$final_Annotations, size_factor = 3)




#Check out the marker expression in the Yalcinbas SCE object, using the final_Annotations

#Custom bubble plot, gets mean expression per cluster for a gene, plots the z-score of that across the clusters
get_bubble_plot_sce = function(sce_object, top_markers, sample_name, group_col = "meta_cluster"){
  # Extract expression data and metadata
  expr_data <- assay(sce_object, "cpm")[top_markers, ]
  metadata <- colData(sce_object)

  # Convert to data frame for plotting
  plot_data <- as.data.frame(t(expr_data)) %>%
    tibble::rownames_to_column("cell_id") %>%
    cbind(group_var = metadata[[group_col]]) %>%
    tidyr::pivot_longer(cols = -c(cell_id, group_var), names_to = "gene", values_to = "expression")

  # Calculate mean expression and percent expressing per cluster
  summary_data <- plot_data %>%
    group_by(gene, group_var) %>%
    summarise(
      mean_expression = mean(expression),
      pct_expressing = sum(expression > 0) / n() * 100,
      .groups = "drop"
    ) %>%
    # Calculate z-score of mean_expression per gene across clusters
    group_by(gene) %>%
    mutate(mean_expression_zscore = scale(mean_expression)[,1]) %>%
    ungroup()

  # Set factor levels to control axis order
  summary_data$gene <- factor(summary_data$gene, levels = top_markers)
  summary_data$group_var <- factor(summary_data$group_var, 
                                      levels = sort(unique(summary_data$group_var)))

  # Create bubble plot
  p1 = ggplot(summary_data, aes(x = gene, y = group_var, size = mean_expression, color = pct_expressing)) +
    geom_point() +
    scale_color_gradient(low = "lightgrey", high = "red") +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() + ggtitle(sample_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = group_col, size = "Mean Expression", color = "% Expressing")

  p2 = ggplot(summary_data, aes(x = gene, y = group_var, size = pct_expressing, color = mean_expression_zscore)) +
    geom_point() +
    scale_color_gradient2(low = "blue", mid = 'white', high = "red", name = "Mean Exp. z-score") +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() + ggtitle(sample_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = group_col, size = "% Expressing", color = "Mean Exp. z-score")
  return(list(p1, p2))
}


#Add CPMs
assay(yalcinbas_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(yalcinbas_sce, 'counts'))

p_bubble = get_bubble_plot_sce(yalcinbas_sce, 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1','TAC3', 'GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1'),
 sample_name = "Yalcinbas Habenula", group_col = "final_Annotations")

p_bubble[[1]]
p_bubble[[2]]


#Compute human markers
markers_human = MetaMarkers::compute_markers(assay(yalcinbas_sce, "cpm"), yalcinbas_sce$final_Annotations)

markers_human %>% group_by(cell_type) %>% slice_max(order_by = auroc, n = 25) %>% View()


