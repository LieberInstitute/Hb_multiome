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
new_data_path = here('processed-data', '98_external_Hb_comparisons', '06_metaMarkers_pandey_2018')
#Path to already generated data
prev_data_path = here('processed-data', '98_external_Hb_comparisons', '05_initial_qc_pandey_2018')
#Path to plot directory
plot_path = here('plots', '98_external_Hb_comparisons', '06_metaMarkers_pandey_2018')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


#Source the bubble plot functions
source(here('code','98_external_Hb_comparisons', 'bubble_plot_functions.R'))


#Load up the full SCE object, contains the metacluster annotations
all_donor_sce = readRDS(paste0(prev_data_path, '/all_donor_sce_with_denovo_clusters.rds'))

#Generate MetaMarkers from the metaclusters
#Split into individual donor SCE objects and compute intial markers for the metaclusters
table(all_donor_sce$study_id)

sce_zebAdult1 <- all_donor_sce[, all_donor_sce$study_id == "Adult1"]
sce_zebAdult2 <- all_donor_sce[, all_donor_sce$study_id == "Adult2"]
#sce_zebLarva <- all_donor_sce[, all_donor_sce$study_id == "Larva"]


markers_adult1 = compute_markers(assay(sce_zebAdult1, "cpm"), sce_zebAdult1$meta_cluster)
markers_adult2 = compute_markers(assay(sce_zebAdult2, "cpm"), sce_zebAdult2$meta_cluster)
#markers_larva = compute_markers(assay(sce_zebLarva, "cpm"), sce_zebLarva$meta_cluster)

head(markers_adult1)


#Save markers
export_markers(markers_adult1, paste0(new_data_path, '/markers_adult1_meta_clusters_markers.csv'))
export_markers(markers_adult2, paste0(new_data_path, '/markers_adult2_meta_clusters_markers.csv'))
#export_markers(markers_larva, paste0(new_data_path, '/markers_larva_meta_clusters_markers.csv'))


#Load up markers and get the metaMarkers 
pandey_zebrafish_hab_markers = list(
    markers_adult1 = read_markers(paste0(new_data_path, '/markers_adult1_meta_clusters_markers.csv.gz')),
    markers_adult2 = read_markers(paste0(new_data_path, '/markers_adult2_meta_clusters_markers.csv.gz'))
    #markers_larva = read_markers(paste0(new_data_path, '/markers_larva_meta_clusters_markers.csv.gz'))
    
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


pandey_zebrafish_hab_metaM %>% group_by(cell_type) %>% slice_min(rank, n = 50) %>% 
  select(cell_type, rank, gene, human_gene_ortholog, recurrence, auroc) %>% 
  View()



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
DimPlot(all_donor_seurat, reduction = "umap") + NoLegend()
#Seurat object without the larval data, just double checking the developmental data is not driving some of the marker trends
#specifically the excitatory/inhibitory population
#no_larva_seurat = subset(all_donor_seurat, subset = study_id != 'Larva')




#In the introduction to the paper, it lists markers for the 3 defined domains
#nptx2a - dorsolateral domain
#gpr151 and pou4f1 dorsomedial domain
#aoc1 ventral domain
zeb_region_markers = c('nptx2a','gpr151','pou4f1','aoc1')
zeb_region_markers = toupper(zeb_region_markers)
p_bubble = get_bubble_plot(all_donor_seurat, zeb_region_markers, 'Zebrafish: Hab region markers')
p_bubble[[1]]
p_bubble[[2]]
#From these markers
#dorsolateral: 3
#dorsomedial: 4, 5, 6
#ventral: 1, 2
# unknown: 7, 8


#Markers used in the original Wallace 2019 paper in Figure 1
custom_markers = c('Tac2', 'Slc17a7', 'Slc17a6', 'Snap25', 'Gap43', 'Slc6a11', 'Cldn5', 'Abcc9', 'Pdgfrb', 'Cx3cr1', 'Mrc1', 'Col3a1', 'Gpr17', 'Mog', 'Olig1', 'Pdgfra')
custom_markers = toupper(custom_markers)
custom_markers = hu_mu_zf_ortholog_df %>% filter(`Gene name` %in% custom_markers) %>% pull(`Zebrafish gene name`) %>% toupper()
custom_markers = custom_markers[custom_markers != '']
p_bubble = get_bubble_plot(all_donor_seurat, custom_markers, 'Zebrafish: Wallace mouse markers')
p_bubble[[1]]
p_bubble[[2]]

#From these markers, high counts of snap25a mask the color scale
#non-neuronal: 8, no snap25 and expresses cldn5a, mrc1a, SLC6A11B, all at low levels
#VGLUT1: 3, 4, 6, though mostly 3 and 4
#VGLUT2: pretty much all the neuronal clusters, 1, 2, 3, 4, 5, 6, 7

#Check the DE of the VGLUT genes, so they're not strong markers, would need a strong outgroup if they all express them at low levels
pandey_zebrafish_hab_metaM %>% filter(gene %in% c('SLC17A6A','SLC17A6B','SLC17A7A','SLC17A7B')) %>% group_by(cell_type) %>% 
  select(cell_type, rank, gene, human_gene_ortholog, recurrence, auroc) %>% 
  View()

#Our human Hab panel
custom_markers = c('Tac3', 'Tac2','Gpr151','Pou4f1','Mbp')
custom_markers = toupper(custom_markers)
custom_markers = hu_mu_zf_ortholog_df %>% filter(`Gene name` %in% custom_markers) %>% pull(`Zebrafish gene name`) %>% toupper()
custom_markers = custom_markers[custom_markers != '']

p_bubble = get_bubble_plot(all_donor_seurat, custom_markers, 'Zebrafish: Human Hab panel markers')
p_bubble[[1]]
p_bubble[[2]]

#Non-neuronal: 8, again with MBPa
#Tac3a: 4, 5 , might be the most aligned with mouse and human medial



#Cholinergic and substance P markers
custom_markers = c('Chat', 'Slc18a3', 'Slc5a7','Tac1', 'Slc17a7', 'Slc17a6')
custom_markers = toupper(custom_markers)
custom_markers = hu_mu_zf_ortholog_df %>% filter(`Gene name` %in% custom_markers) %>% pull(`Zebrafish gene name`) %>% toupper()
custom_markers = custom_markers[custom_markers != '']
p_bubble = get_bubble_plot(all_donor_seurat, custom_markers, 'Zebrafish: Cholinergic and substance P markers')
p_bubble[[1]]
p_bubble[[2]]

#Cholinergic, so there's no expression of chata and chatb is not in the gene annotations
#slc5a7a: 4, 5, only cholinergic marker that seems to have expression
#supstanc P: 3



#Adult Zebrafish markers used in Pandey Figures
zeb_custom_adult_markers = c('tac3a','adrb2a','gng2','cbln2b','trh','lrrtm1','wnt7aa','adcyap1a','igf2a',
'pvalb7','sox1b','tubb5','gad2','cntnap2a', 'rgs5b','cd82a','zgc:173443','her4.3')
zeb_custom_adult_markers = toupper(zeb_custom_adult_markers)
p_bubble = get_bubble_plot(all_donor_seurat, zeb_custom_adult_markers, 'Zebrafish all 3 samples')
p_bubble[[1]]
p_bubble[[2]]

zeb_custom_adult_markers = c('tac3a','tac3b','adrb2a','adrb2b','adcyap1a','adcyap1b',
'igf2a')
zeb_custom_adult_markers = toupper(zeb_custom_adult_markers)
p_bubble = get_bubble_plot(all_donor_seurat, zeb_custom_adult_markers, 'Zebrafish all 3 samples')
p_bubble[[1]]
p_bubble[[2]]

#The right hab marker used was Tac3a, corresponds to 4 and 5
#The left hab markers used were ADYAP1A and IGF2A, correspond to cluster 3

#Larval zebrafish markers used in Pandey figures
#zeb_custom_larva_markers = c('murcb','adrb2a','spx','cbln2b','c1ql4b','lrrtm1','pcdh7b','wnt7aa','adcyap1a', 'igf2a',
#'ppp1r1c','sox1a','htr1aa','tubb5','gad2','kiss1', 'epcam')
#zeb_custom_larva_markers = toupper(zeb_custom_larva_markers)
#p_bubble = get_bubble_plot(all_donor_seurat, zeb_custom_larva_markers, 'Zebrafish all 3 samples')
#p_bubble[[1]]
#p_bubble[[2]]



#Top 10 metamarkers per metacluster
custom_meta_markers = pandey_zebrafish_hab_metaM %>% group_by(cell_type) %>% slice_min(rank, n = 10) %>% pull(gene)
custom_meta_markers = unique(custom_meta_markers)
p_bubble = get_bubble_plot(all_donor_seurat, custom_meta_markers, 'Zebrafish: Top 10 markers')
p_bubble[[1]]
p_bubble[[2]]



custom_markers = c('gad1','gad2','slc32a1','slc17a6','slc17a7', 'OPRM1', 'OPRD1','OPRK1')
custom_markers = toupper(custom_markers)
custom_markers = hu_mu_zf_ortholog_df %>% filter(`Gene name` %in% custom_markers) %>% pull(`Zebrafish gene name`) %>% toupper()
custom_markers = custom_markers[custom_markers != '']
p_bubble = get_bubble_plot(all_donor_seurat, custom_markers, 'Zebrafish: Inhibitory/Excitatory and opioid markers')
p_bubble[[1]]
p_bubble[[2]]



#Inhibitory: 7, expresses GAD1, GAD2, and VGAT, low levels of VGLUT2



#Zebrafish specific cholinergic markers, found from https://www-nature-com.proxy1.library.jhu.edu/articles/s41598-020-72524-3
#Presynaptic markers: chata, chatb, hacta, hactb, vachta, vachtb, ache
#nicotinic ach reseptors: chrna3, chrna7, chrna2a
#cholinergic channels?: slc5a7
#supstance p markers: tac1, tacr1a, tacr1b

custom_chol_markers = c('chata', 'chatb', 'hacta', 'hactb', 'vachta', 'vachtb', 'ache', 'chrna3', 'chrna7', 'chrna2a',
'slc5a7a', 'slc5a7','tac1', 'tacr1a', 'tacr1b')
custom_chol_markers = toupper(custom_chol_markers)
p_bubble = get_bubble_plot(all_donor_seurat, custom_chol_markers, 'Zebrafish: Cholinergic and substance P markers')
p_bubble[[1]]
p_bubble[[2]]



pandey_zebrafish_hab_metaM %>% filter(gene == 'SLC6A1B') %>% 
  View()


#Taking all of the above together, here are initial annotations
#1: ventral
#2: ventral, top markers seem to include immediate early genes
#3: dorsolateral_TAC1_left, BDNF is also a marker
#4: dorsomedial_TAC3A_cholinergic_right
#5: dorsomedial_TAC3A_cholinergic_right, stronger expression of tac3 and slc5A7 gene, also SLC6A1 (GAT1) is a marker
#6: dorsomedial_neuron
#7: inhibitory_gap43, VGAT, GAD1, and GAD2 are all markers
#8: non-neuronal, alot of ribosomal genes are top markers, might be general mush


#metacluster annotations from all the above
meta_annot_vec = c( 'ventral' = 'meta_cluster1',
                    'ventral_immediate_early' = 'meta_cluster2',
                    'dorsolateral_left_subP_BDNF' = 'meta_cluster3',
                    'dorsomedial_right_cholinergic' = 'meta_cluster4',
                    'dorsomedial_right_cholinergic_GAT1' = 'meta_cluster5',
                    'dorsomedial_neuron' = 'meta_cluster6',
                    'inhibitory_gap43' = 'meta_cluster7',
                    'non_neuronal' = 'meta_cluster8',
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

p6 = DimPlot(all_donor_seurat , reduction = "umap", group.by = 'study_id', label = TRUE) + 
  ggtitle('Pandey 2018: sample batch')
p6




#With the annotated clusters, show and save the bubble plot with the same genes as the mouse data


custom_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1', 'TAC2','GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1')
custom_markers = toupper(custom_markers)
custom_markers = hu_mu_zf_ortholog_df %>% filter(`Gene name` %in% custom_markers) %>% pull(`Zebrafish gene name`) %>% toupper()
custom_markers = custom_markers[custom_markers != '']
custom_markers

p_bubble = get_bubble_plot(all_donor_seurat, 
  top_markers = custom_markers, sample_name = "Zebrafish Habenula", group_col = "meta_clust_celltype_annot")
p_bubble[[1]]
p_bubble[[2]]

pdf(paste0(plot_path, '/Zebrafish_Hab_marker_bubbles_meanExp.pdf'), width = 10, height = 8)
p_bubble[[1]]
dev.off()
pdf(paste0(plot_path, '/Zebrafish_Hab_marker_bubbles_Zscore_meanExp.pdf'), width = 10, height = 8)
p_bubble[[2]]
dev.off()



#And the regional zebrafish markers
zeb_region_markers = c('nptx2a','gpr151','pou4f1','aoc1')
zeb_region_markers = toupper(zeb_region_markers)
p_bubble = get_bubble_plot(all_donor_seurat, 
  top_markers = zeb_region_markers, 
  sample_name = "Zebrafish Habenula", group_col = "meta_clust_celltype_annot")
p_bubble[[1]]
p_bubble[[2]]

pdf(paste0(plot_path, '/Zebrafish_Hab_regional_marker_bubbles_meanExp.pdf'), width = 10, height = 8)
p_bubble[[1]]
dev.off()
pdf(paste0(plot_path, '/Zebrafish_Hab_regional_marker_bubbles_Zscore_meanExp.pdf'), width = 10, height = 8)
p_bubble[[2]]
dev.off()


##########################
#Look for potential lateral habenula GABAergic cells
#########################


custom_markers = c('GAD1A', 'GAD1B',  'GAD2',  'SLC32A1', 'SLC17A6A', 'SLC17A6B', 'SLC17A7A', 'SLC17A7B')


p_bubble_gaba_glut_markers = get_bubble_plot(all_donor_seurat , custom_markers, 'Pendey Zebrafish', 
group_col = 'meta_clust_celltype_annot')
p_bubble_gaba_glut_markers


#Here, the inhibitory_gap43 population clearly has the most distince GAD expression
zeb_gad1a_p = FeaturePlot(all_donor_seurat, features = "GAD1A", reduction = "umap", pt.size = 1, slot = 'scale.data') +
  scale_color_gradient(low = "white", high = "red", name = 'z-score')

zeb_gad1b_p = FeaturePlot(all_donor_seurat, features = "GAD1B", reduction = "umap", pt.size = 1, slot = 'scale.data') +
  scale_color_gradient(low = "white", high = "red", name = 'z-score')

zeb_gad2_p = FeaturePlot(all_donor_seurat, features = "GAD2", reduction = "umap", pt.size = 1, slot = 'scale.data') +
  scale_color_gradient(low = "white", high = "red", name = 'z-score')

zeb_vgat_p = FeaturePlot(all_donor_seurat, features = "SLC32A1", reduction = "umap", pt.size = 1, slot = 'scale.data') +
  scale_color_gradient(low = "white", high = "red", name = 'z-score')



#Look at the co-expression of specific genes
# Get expression data
umap_data <- as.data.frame(Embeddings(all_donor_seurat, reduction = "umap"))
umap_data$GAD2 <- FetchData(all_donor_seurat, vars = "GAD2", slot = "data")[, 1]
umap_data$SLC17A6A <- FetchData(all_donor_seurat, vars = "SLC17A6A", slot = "data")[, 1]

# Create a coexpression category
umap_data$coexpression <- ifelse(umap_data$GAD2 > 0 & umap_data$SLC17A6A > 0, "Both",
                                  ifelse(umap_data$GAD2 > 0, "GAD2 only",
                                         ifelse(umap_data$SLC17A6A > 0, "SLC17A6A only", "Neither")))

gad2_vglut2_p = ggplot(umap_data, aes(x = umap_1, y = umap_2, color = coexpression)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Both" = "purple", "GAD2 only" = "red", "SLC17A6A only" = "blue", "Neither" = "lightgrey")) +
  theme_bw() +
  labs(title = "GAD2 and SLC17A6A Co-expression: Pendey Zebrafish")


p_bubble_gaba_glut_markers[[1]]
p_bubble_gaba_glut_markers[[2]]

zeb_gad1a_p 
zeb_gad1b_p
zeb_gad2_p
zeb_vgat_p

gad2_vglut2_p 

ggsave(p_bubble_gaba_glut_markers[[1]], filename = 'pendey_zeb_gaba_glut_meta_annots_bubble.pdf', path = plot_path,
device = 'pdf', width = 10, height = 8)

ggsave(p_bubble_gaba_glut_markers[[2]], filename = 'pendey_zeb_gaba_glut_zscore_meta_annots_bubble.pdf', path = plot_path,
device = 'pdf', width = 10, height = 8)

ggsave(zeb_gad1a_p , filename = 'pendey_zeb_gad1a_exp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(zeb_gad1b_p , filename = 'pendey_zeb_gad1b_exp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(zeb_gad2_p , filename = 'pendey_zeb_gad2_exp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(zeb_vgat_p , filename = 'pendey_zeb_vgat_exp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

ggsave(gad2_vglut2_p , filename = 'pendey_zeb_gad2_vglut2_coexp_umap.pdf', path = plot_path,
device = 'pdf', width = 6, height = 4)

#Coexpression plots of the GABA-Glut genes

all_donor_sce$meta_clust_celltype_annot = unname(meta_annot_vec[all_donor_sce$meta_cluster])
logcounts(all_donor_sce) = log1p(assay(all_donor_sce, "cpm"))
#Plotting co-expression of excitatory and inhibitory markers

pairwise_coexpression <- function(mat, genes, cluster_name) {
  # mat: genes x cells matrix for one cluster
  
  detected <- mat[genes, , drop = FALSE] > 0
  
  res <- expand.grid(gene1 = genes, gene2 = genes, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(percent = mean(detected[gene1, ] & detected[gene2, ]) * 100) %>%
    ungroup() %>%
    mutate(cluster = cluster_name)
  
  res
}


#Full dataset. Of note, this does not use the corrected counts, the corrected counts were not saved with this version of the data
genes = c('GAD1A', 'GAD1B',  'GAD2',  'SLC32A1', 'SLC17A6A', 'SLC17A6B', 'SLC17A7A', 'SLC17A7B')

#genes <- c('Slc32a1',"Gad1","Gad2","Slc17a6", "Slc17a7")
expr_mat <- assay(all_donor_sce, "logcounts")

cluster_to_annotate = "meta_clust_celltype_annot"
clusters <- unique(colData(all_donor_sce)[[cluster_to_annotate]])

coexp_df <- lapply(clusters, function(cl) {
  cells <- colData(all_donor_sce)[[cluster_to_annotate]] == cl
  mat_sub <- expr_mat[, cells, drop = FALSE]
  pairwise_coexpression(mat_sub, genes, cluster_name = cl)
}) %>%
  bind_rows()

coexp_df$gene1 <- factor(coexp_df$gene1, levels = genes)
coexp_df$gene2 <- factor(coexp_df$gene2, levels = rev(genes))


#Edited the original code to have grey be between 0-5%, previously the very low percentages were difficult to see.
p <- ggplot(coexp_df, aes(x = gene1, y = gene2, fill = percent)) +
  geom_tile(color = "grey70", linewidth = 0.3) +
  facet_wrap(~ cluster, nrow = 5) +
  scale_fill_gradientn(
    colours = c("grey95","#f1e2c6", "#f1e2c6", "#f0c94a", "#df8b27", "#d92523", "#8b0d19"),
    values = c(0, 0.05, 0.20, 0.40, 0.60, 0.8, 1),
    limits = c(0, 100),
    breaks = c( 5, 20, 40, 60, 80, 100),
    name = "Percent of cells expressing\ntwo genes"
  ) +
  coord_equal() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  xlab(NULL) +
  ylab(NULL) + ggtitle("Co-expression of GABA and Glut markers in Zebrafish Pendey dataset")

print(p)

ggsave(path = plot_path, filename = 'GABA_Glut_coexpression_Pendey_zeb.pdf', plot = p, 
device = 'pdf', width = 14, height = 10, useDingbats = FALSE)



























