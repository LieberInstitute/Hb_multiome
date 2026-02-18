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

path_to_orthologs = here('processed-data', '98_external_Hb_comparisons', 'human_mouse_zebrafish_orthologs.txt.gz')

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


#Add the human orthologs to the zebrafish metamarkers
#`Gene name` column is the human gene name
hu_mu_zf_ortholog_df = data.table::fread(path_to_orthologs)
table(hu_mu_zf_ortholog_df$`Mouse homology type`)
#Get the 1to1 orthologs across the three species
hu_mu_zf_ortholog_df <- hu_mu_zf_ortholog_df %>% filter(`Mouse homology type` == 'ortholog_one2one' )
hu_mu_zf_ortholog_df  = hu_mu_zf_ortholog_df %>% filter(!duplicated(`Gene name`))

dim(hu_mu_zf_ortholog_df )
View(hu_mu_zf_ortholog_df)
#The zebrafish gene names are all uppercase in the data, need to match the lowercase for the orthologtable
head(pandey_zebrafish_hab_metaM )
pandey_zebrafish_hab_metaM$lowercase_zeb_gene = tolower(pandey_zebrafish_hab_metaM$gene)

#Add the human ortholog
index = match(pandey_zebrafish_hab_metaM$lowercase_zeb_gene, hu_mu_zf_ortholog_df$`Zebrafish gene name` )
pandey_zebrafish_hab_metaM$human_gene_ortholog = hu_mu_zf_ortholog_df$`Gene name`[index]


pandey_zebrafish_hab_metaM %>% group_by(cell_type) %>% slice_min(rank, n = 20) %>% 
  select(cell_type, rank, gene, human_gene_ortholog, recurrence, auroc) %>% 
  View()



#Custom bubble plot function

get_bubble_plot = function(seurat_object, top_markers, sample_name, group_col = 'meta_cluster'){
  # Extract expression data and metadata
  expr_data <- FetchData(seurat_object, vars = top_markers, slot = "data")
  metadata <- seurat_object@meta.data

  # Combine into a data frame
  plot_data <- cbind(expr_data, meta_cluster = metadata[[group_col]]) %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = -meta_cluster, names_to = "gene", values_to = "expression")

  # Calculate mean expression and percent expressing per cluster
  summary_data <- plot_data %>%
    group_by(gene, meta_cluster) %>%
    summarise(
      mean_expression = mean(expression),
      pct_expressing = sum(expression > 0) / n() * 100,
      .groups = "drop"
    )

  # Set factor levels to control axis order
  summary_data$gene <- factor(summary_data$gene, levels = top_markers)
  summary_data$meta_cluster <- factor(summary_data$meta_cluster, 
                                      levels = sort(unique(summary_data$meta_cluster)))

  # Create bubble plot
  p1 = ggplot(summary_data, aes(x = gene, y = meta_cluster, size = mean_expression, color = pct_expressing)) +
    geom_point() +
    scale_color_gradient(low = "lightgrey", high = "red") +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() + ggtitle(sample_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = "Meta Cluster", size = "Mean Expression", color = "% Expressing")
  return(p1)
}


#Switch the full dataset to Seurat for the bubble plots
all_donor_seurat = as.Seurat(all_donor_sce, counts = "counts", data = "cpm")


#Markers used in the original Wallace 2019 paper in Figure 1
custom_markers = c('Tac2', 'Slc17a7', 'Slc17a6', 'Snap25', 'Gap43', 'Slc6a11', 'Cldn5', 'Abcc9', 'Pdgfrb', 'Cx3cr1', 'Mrc1', 'Col3a1', 'Gpr17', 'Mog', 'Olig1', 'Pdgfra')
custom_markers = toupper(custom_markers)
custom_markers = hu_mu_zf_ortholog_df %>% filter(`Gene name` %in% custom_markers) %>% pull(`Zebrafish gene name`) %>% toupper()
custom_markers = custom_markers[custom_markers != '']
p_bubble = get_bubble_plot(all_donor_seurat, custom_markers, 'Zebrafish all 3 samples')
p_bubble



custom_markers = c('Tac3', 'Tac2','Gpr151','Pou4f1','Mbp')
custom_markers = toupper(custom_markers)
custom_markers = hu_mu_zf_ortholog_df %>% filter(`Gene name` %in% custom_markers) %>% pull(`Zebrafish gene name`) %>% toupper()
custom_markers = custom_markers[custom_markers != '']

p_bubble = get_bubble_plot(all_donor_seurat, custom_markers, 'Zebrafish all 3 samples')
p_bubble


custom_markers = c('Chat', 'Slc18a3', 'Slc5a7','Tac1', 'Slc17a7', 'Slc17a6')
custom_markers = toupper(custom_markers)
custom_markers = hu_mu_zf_ortholog_df %>% filter(`Gene name` %in% custom_markers) %>% pull(`Zebrafish gene name`) %>% toupper()
custom_markers = custom_markers[custom_markers != '']
p_bubble = get_bubble_plot(all_donor_seurat, custom_markers, 'Zebrafish all 3 samples')
p_bubble




























