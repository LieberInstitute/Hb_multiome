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

#Get the 1to1 orthologs across human and mouse, will likely need to relax for zebrafish, I think it's still reasonable to work with the many-to-many orthologs there
hu_mu_zf_ortholog_df <- hu_mu_zf_ortholog_df %>% filter(`Mouse homology type` == 'ortholog_one2one' )

#To keep a record of the many-to-many zebrafish orthologs, filter on the unique combo of all three gene name columns
hu_mu_zf_ortholog_df$comb_species_gene = paste(paste(hu_mu_zf_ortholog_df$`Gene name`, hu_mu_zf_ortholog_df$`Mouse gene name`, sep = '_'), hu_mu_zf_ortholog_df$`Zebrafish gene name`, sep = '_' )
hu_mu_zf_ortholog_df  = hu_mu_zf_ortholog_df %>% filter(!duplicated(comb_species_gene))

dim(hu_mu_zf_ortholog_df )
hu_mu_zf_ortholog_df %>% 
  select(`Gene name`, `Mouse gene name`, `Mouse homology type`, `Zebrafish gene name`, `Zebrafish homology type`) %>%
  View()

#The zebrafish gene names are all uppercase in the data, need to match the lowercase for the orthologtable
head(pandey_zebrafish_hab_metaM )
pandey_zebrafish_hab_metaM$lowercase_zeb_gene = tolower(pandey_zebrafish_hab_metaM$gene)

#Add the human ortholog
index = match(pandey_zebrafish_hab_metaM$lowercase_zeb_gene, hu_mu_zf_ortholog_df$`Zebrafish gene name` )
pandey_zebrafish_hab_metaM$human_gene_ortholog = hu_mu_zf_ortholog_df$`Gene name`[index]


pandey_zebrafish_hab_metaM %>% group_by(cell_type) %>% slice_min(rank, n = 20) %>% 
  select(cell_type, rank, gene, human_gene_ortholog, recurrence, auroc) %>% 
  View()


#Custom bubble plot, gets mean expression per cluster for a gene, plots the z-score of that across the clusters
get_bubble_plot = function(seurat_object, top_markers, sample_name, group_col = "meta_cluster"){
  # Extract expression data and metadata
  expr_data <- FetchData(seurat_object, vars = top_markers, layer = "data")
  metadata <- seurat_object@meta.data

  # Combine into a data frame
  plot_data <- cbind(expr_data, group_var = metadata[[group_col]]) %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols = -group_var, names_to = "gene", values_to = "expression")

  # Calculate mean expression and percent expressing per cluster
  summary_data <- plot_data %>% filter(group_var != 'outliers') %>%
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
  p1 = ggplot(summary_data, aes(x = gene, y = group_var, size = pct_expressing, color = mean_expression)) +
    geom_point() +
    scale_color_gradient2(low = "white", high = "red", name = "Mean Expression") +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() + ggtitle(sample_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = group_col, size = "% Expressing", color = "Mean Expression")
    
  p2 = ggplot(summary_data, aes(x = gene, y = group_var, size = pct_expressing, color = mean_expression_zscore)) +
    geom_point() +
    scale_color_gradient2(low = "blue", mid = 'white', high = "red", name = "Mean Exp. z-score") +
    scale_size_continuous(range = c(2, 8)) +
    theme_minimal() + ggtitle(sample_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Gene", y = group_col, size = "% Expressing", color = "Mean Exp. z-score")
  return(list(p1, p2))
}


#Switch the full dataset to Seurat for the bubble plots
all_donor_seurat = as.Seurat(all_donor_sce, counts = "counts", data = "cpm")
#HVGs
all_donor_seurat <- FindVariableFeatures(all_donor_seurat, selection.method = "vst", nfeatures = 2000)

#Scale data
all.genes <- rownames(all_donor_seurat)
all_donor_seurat <- ScaleData(all_donor_seurat, features = all.genes)

#PCA
all_donor_seurat  <- RunPCA(all_donor_seurat , features = VariableFeatures(object = all_donor_seurat ))
DimPlot(all_donor_seurat, reduction = "pca") + NoLegend()

#UMAP
all_donor_seurat  <- RunUMAP(all_donor_seurat , dims = 1:20)


#In the introduction to the paper, it lists markers for the 3 defined domains
#nptx2a - dorsolateral domain
#gpr151 and pou4f1 dorsomedial domain
#aoc1 ventral domain
zeb_region_markers = c('nptx2a','gpr151','pou4f1','aoc1')
zeb_region_markers = toupper(zeb_region_markers)
p_bubble = get_bubble_plot(all_donor_seurat, zeb_region_markers, 'Zebrafish all 3 samples')
p_bubble[[1]]
p_bubble[[2]]
#From these markers
#dorsolateral: 8
#dorsomedial: 1, 2, 3, 4, 7, 9
#ventral: 5, 6
# unknown: 10 


#Markers used in the original Wallace 2019 paper in Figure 1
custom_markers = c('Tac2', 'Slc17a7', 'Slc17a6', 'Snap25', 'Gap43', 'Slc6a11', 'Cldn5', 'Abcc9', 'Pdgfrb', 'Cx3cr1', 'Mrc1', 'Col3a1', 'Gpr17', 'Mog', 'Olig1', 'Pdgfra')
custom_markers = toupper(custom_markers)
custom_markers = hu_mu_zf_ortholog_df %>% filter(`Gene name` %in% custom_markers) %>% pull(`Zebrafish gene name`) %>% toupper()
custom_markers = custom_markers[custom_markers != '']
p_bubble = get_bubble_plot(all_donor_seurat, custom_markers, 'Zebrafish all 3 samples')
p_bubble[[1]]
p_bubble[[2]]

#From these markers, high counts of snap25a mask the color scale
#non-neuronal: 10, no snap25 and expresses cldn5a, mrc1a, abss9, all at low levels
#VGLUT1: 7, 8
#VGLUT2: 1, 2, 3, 7, 8

#Our human Hab panel
custom_markers = c('Tac3', 'Tac2','Gpr151','Pou4f1','Mbp')
custom_markers = toupper(custom_markers)
custom_markers = hu_mu_zf_ortholog_df %>% filter(`Gene name` %in% custom_markers) %>% pull(`Zebrafish gene name`) %>% toupper()
custom_markers = custom_markers[custom_markers != '']

p_bubble = get_bubble_plot(all_donor_seurat, custom_markers, 'Zebrafish all 3 samples')
p_bubble[[1]]
p_bubble[[2]]

#Non-neuronal: 10, again with MBPa
#Tac3a: 1, 2, might be the most aligned with mouse and human medial



#Cholinergic and substance P markers
custom_markers = c('Chat', 'Slc18a3', 'Slc5a7','Tac1', 'Slc17a7', 'Slc17a6')
custom_markers = toupper(custom_markers)
custom_markers = hu_mu_zf_ortholog_df %>% filter(`Gene name` %in% custom_markers) %>% pull(`Zebrafish gene name`) %>% toupper()
custom_markers = custom_markers[custom_markers != '']
p_bubble = get_bubble_plot(all_donor_seurat, custom_markers, 'Zebrafish all 3 samples')
p_bubble[[1]]
p_bubble[[2]]

#Cholinergic, so there's no expression of chata and chatb is not in the gene annotations
#slc5a7a: 1, 2, only cholinergic marker that seems to have expression
#supstanc P: 4, 8



#Adult Zebrafish markers used in Pandey Figures
zeb_custom_adult_markers = c('tac3a','adrb2a','gng2','cbln2b','trh','lrrtm1','wnt7aa','adcyap1a',
'pvalb7','sox1b','tubb5','gad2','cntnap2a', 'rgs5b','cd82a','zgc:173443','her4.3')
zeb_custom_adult_markers = toupper(zeb_custom_adult_markers)
p_bubble = get_bubble_plot(all_donor_seurat, zeb_custom_adult_markers, 'Zebrafish all 3 samples')
p_bubble[[1]]
p_bubble[[2]]

#Larval zebrafish markers used in Pandey figures
zeb_custom_larva_markers = c('murcb','adrb2a','spx','cbln2b','c1ql4b','lrrtm1','pcdh7b','wnt7aa','adcyap1a',
'ppp1r1c','sox1a','htr1aa','tubb5','gad2','kiss1', 'epcam')
zeb_custom_larva_markers = toupper(zeb_custom_larva_markers)
p_bubble = get_bubble_plot(all_donor_seurat, zeb_custom_larva_markers, 'Zebrafish all 3 samples')
p_bubble[[1]]
p_bubble[[2]]



#Top 10 metamarkers per metacluster
custom_meta_markers = pandey_zebrafish_hab_metaM %>% group_by(cell_type) %>% slice_min(rank, n = 10) %>% pull(gene)
custom_meta_markers = unique(custom_meta_markers)
p_bubble = get_bubble_plot(all_donor_seurat, custom_meta_markers, 'Zebrafish all 3 samples')
p_bubble[[1]]
p_bubble[[2]]
#Only real thing that stands out is inhibitory markers in cluster 3
#This population expresses gad1, gad2, vglut1, vglut2 


custom_markers = c('gad1','gad2','slc32a1','slc17a6','slc17a7', 'OPRM1', 'OPRD1','OPRK1')
custom_markers = toupper(custom_markers)
custom_markers = hu_mu_zf_ortholog_df %>% filter(`Gene name` %in% custom_markers) %>% pull(`Zebrafish gene name`) %>% toupper()
custom_markers = custom_markers[custom_markers != '']
p_bubble = get_bubble_plot(all_donor_seurat, custom_markers, 'Zebrafish all 3 samples')
p_bubble[[1]]
p_bubble[[2]]

#Inhibitory: 3, expresses GAD1, GAD2, and VGAT, low levels of VGLUT2



#Zebrafish specific cholinergic markers, found from https://www-nature-com.proxy1.library.jhu.edu/articles/s41598-020-72524-3
#Presynaptic markers: chata, chatb, hacta, hactb, vachta, vachtb, ache
#nicotinic ach reseptors: chrna3, chrna7, chrna2a
#cholinergic channels?: slc5a7
#supstance p markers: tac1, tacr1a, tacr1b

custom_chol_markers = c('chata', 'chatb', 'hacta', 'hactb', 'vachta', 'vachtb', 'ache', 'chrna3', 'chrna7', 'chrna2a',
'slc5a7a', 'slc5a7','tac1', 'tacr1a', 'tacr1b')
custom_chol_markers = toupper(custom_chol_markers)
p_bubble = get_bubble_plot(all_donor_seurat, custom_chol_markers, 'Zebrafish all 3 samples')
p_bubble[[1]]
p_bubble[[2]]


#Taking all of the above together, here are initial annotations
#1: dorsomedial_cholinergic_tac3a_1
#2: dorsomedial_cholinergic_tac3a_2
#3: dorsomedial_inhibitory_gap43_vglut2
#4: dorsomedial_subP
#5: ventral_1
#6: ventral_2
#7: dorsomedial_vglut1_vglut2
#8: dorsolateral_subP_vglut1_vglut2
#9: dorsomedial_vglut2
#10: non_neuronal


#metacluster annotations from all the above
meta_annot_vec = c('dorsomedial_cholinergic_tac3a_1' = 'meta_cluster1',
                    'dorsomedial_cholinergic_tac3a_2' = 'meta_cluster2',
                    'dorsomedial_inhibitory_gap43_vglut2' = 'meta_cluster3',
                    'dorsomedial_subP' = 'meta_cluster4',
                    'ventral_1' = 'meta_cluster5',
                    'ventral_2' = 'meta_cluster6',
                    'dorsomedial_vglut1_vglut2' = 'meta_cluster7',
                    'dorsolateral_subP_vglut1_vglut2' = 'meta_cluster8',
                    'dorsomedial_vglut2' = 'meta_cluster9',
                    'non_neuronal' = 'meta_cluster10',
                    'outliers' = 'outliers'  
)

meta_annot_vec  = setNames(names(meta_annot_vec), meta_annot_vec)


all_donor_seurat$meta_clust_celltype_annot = unname(meta_annot_vec[all_donor_seurat$meta_cluster])
table(all_donor_seurat$meta_clust_celltype_annot, all_donor_seurat$meta_cluster)

#Save the metadata as a data.frame to add to the seurat data object later
full_seurat_metadata = all_donor_seurat@meta.data
saveRDS(full_seurat_metadata, paste0(new_data_path, '/pandey_zebrafish_metaclust_celltype_annot_metadata.rds'))



#and a final umap with the annotated metaclusters
p5 = DimPlot(all_donor_seurat , reduction = "umap", group.by = 'meta_clust_celltype_annot', label = TRUE) + 
  ggtitle('Pandey 2018: metacluster celltype annotations')
p5
ggsave(p5, filename = 'pandey_zebrafish_hab_metaCluster_celltype_annot_umap.pdf', path = plot_path,
device = 'pdf', width = 8, height = 7)


































