
#There is some signal for inhibitory neurons within the lateral habenula or coGaba-Glut releasing neurons, currently within LHb.4
#Dig into this a bit more. Check whether these are real co-releasers or could be mislabeled thalamus cells.

library(SingleCellExperiment)
library(Seurat)
library(MetaMarkers)
library(MetaNeighbor)
library(dplyr)
library(ggplot2)
library(scater)
library(scran)
library(here)

here::here()

#Path for new data generated
new_data_path = here('processed-data', '99_donor_cluster_replicability','02_LHb4_investigation')
#Path to plot directory
plot_path = here('plots', '99_donor_cluster_replicability','02_LHb4_investigation')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)

#Source the bubble plot functions
source(here('code','98_external_Hb_comparisons', 'bubble_plot_functions.R'))


multiome_path = here('processed-data', '08_spatial_registration_vs_multiome_snRNA-seq','mid')

#Multiome human data
multiome_sce = readRDS(paste0(multiome_path, '/seurat.norm_counts_CRr_WNN_rnaHarm_atacHarm_k30_C.leiden_lsi_r2_renamed_visium_v5.rds'))
assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))

multiome_sce


#Set up color scale
color_palette_1 = c('Astrocyte' = 'grey','Endo' = 'firebrick','Microglia' = 'darkred',
'Oligo' = 'darkgoldenrod','OPC' = 'cornsilk3')

color_palette_2 = c('Inhib.Thal' = 'dodgerblue','Excit.Thal' = 'indianred2','Thal' = 'lightsteelblue')

color_palette_3 = MetBrewer::met.brewer("Redon", n = 10)
names(color_palette_3) = c('LHb.1','LHb.1.3','LHb.1.3.4','LHb.2.7','LHb.4','LHb.7','MHb.1','MHb.1.2','MHb.2','MHb.3')

color_palette = c(color_palette_1, color_palette_2, color_palette_3)


#Marker plot looking at markers for the basic expected types in the habenula
p_bubble = get_bubble_plot_sce(multiome_sce, 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1','TAC3', 'GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1'),
 sample_name = "Multiome Habenula", group_col = "mid_cluster")

p_bubble[[1]]
p_bubble[[2]]

pdf(paste0(plot_path, '/full_multiome_bubble_Hab_marker_mean.pdf'), width = 10, height = 8)
p_bubble[[1]]
dev.off()

pdf(paste0(plot_path, '/full_multiome_bubble_Hab_marker_zscore_mean.pdf'), width = 10, height = 8)
p_bubble[[2]]
dev.off()

#Very small signal for GAD1, GAD2, VGAT expression in LHb4 and LHb7, but LHb7 is that donor specific cluster


#Initial check on the co-expression of GABA and Glut markers across clusters
#Split by donor too

#Adapted from https://github.com/LieberInstitute/spatial_LS/blob/c68566780bc9a4e1aa3fe3f1c238d424bd4d8447/code/14_annotating_chromium_a-p/02_annotating.R#L150-L197

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
genes <- c('SLC32A1',"GAD1","GAD2","SLC17A6", "SLC17A7")
expr_mat <- assay(multiome_sce, "logcounts")

cluster_to_annotate = "mid_cluster"
clusters <- unique(colData(multiome_sce)[[cluster_to_annotate]])

coexp_df <- lapply(clusters, function(cl) {
  cells <- colData(multiome_sce)[[cluster_to_annotate]] == cl
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
    breaks = c(0, 5, 20, 40, 60, 80, 100),
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
  ylab(NULL) + ggtitle("Co-expression of GABA and Glut markers in full multiome dataset")

print(p)

ggsave(path = plot_path, filename = 'GABA_Glut_coexpression_full_multiome.pdf', plot = p, 
device = 'pdf', width = 12, height = 10, useDingbats = FALSE)



#Check out the signal by donors

donors <- unique(colData(multiome_sce)[['orig.ident']])

for(donor in donors){
  sce_to_plot = multiome_sce[ , multiome_sce$orig.ident == donor]

  genes <- c('SLC32A1',"GAD1","GAD2","SLC17A6", "SLC17A7")
  expr_mat <- assay(sce_to_plot, "logcounts")

  cluster_to_annotate = "mid_cluster"
  clusters <- unique(colData(sce_to_plot)[[cluster_to_annotate]])

  coexp_df <- lapply(clusters, function(cl) {
    cells <- colData(sce_to_plot)[[cluster_to_annotate]] == cl
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
      breaks = c(0, 5, 20, 40, 60, 80, 100),
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
    ylab(NULL) + ggtitle(paste("Co-expression of GABA and Glut markers in", donor))

  print(p)
  
}


#Quick check on the co-expression of CHAT, SLC5A7, and TAC1 for the medial habenula clusters
genes <- c('CHAT',"SLC5A7","TAC1")
expr_mat <- assay(multiome_sce, "logcounts")

cluster_to_annotate = "mid_cluster"
clusters <- unique(colData(multiome_sce)[[cluster_to_annotate]])

coexp_df <- lapply(clusters, function(cl) {
  cells <- colData(multiome_sce)[[cluster_to_annotate]] == cl
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
    breaks = c(0, 5, 20, 40, 60, 80, 100),
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
  ylab(NULL) + ggtitle("Co-expression of Cholinergic and SubP markers in full multiome dataset")

print(p)

ggsave(path = plot_path, filename = 'Chol_subP_coexpression_full_multiome.pdf', plot = p, 
device = 'pdf', width = 12, height = 10, useDingbats = FALSE)



############################
#While the LHb.4 cluster does have some co-expression of GABA and Glut markers, they are super small percentages.
#I'm curious if these are just small subsets of thalamus cells that are present within the LHb.4 cluster
#Two pieces of evidence, the LHb.4 cluster is spatially the one closest to the thalamus
#And the GABA glut co-expression is the strongest in the thalamus clusters.

#So even a small amount of thalamus in the LHb.4 could account for the low levels of co GABA glut

##############################

#For the subset of GABA-Glut cells in a cluster, compare the expression of thalamus markers for those cells

#First though, let's look at the expression of GADs and VGLUTs in just the LHb4 and LHb7 clusters

Lhab_sce = multiome_sce[, multiome_sce$mid_cluster %in% c('LHb.4', 'LHb.7')]

# Get unique donor IDs
donors <- unique(Lhab_sce$orig.ident)

# Create a named list of SCE objects, one per donor
Lhab_sce_list <- lapply(donors, function(donor) {
  Lhab_sce[, Lhab_sce$orig.ident == donor]
})
names(Lhab_sce_list) <- donors

#Convert first to seurats
multiome_seurat_list <- lapply(Lhab_sce_list, function(sce) {
  as.Seurat(sce, counts = "counts", data = "cpm")
})


#Default seurat integration for cross-donors

# Find integration anchors
integration_anchors <- FindIntegrationAnchors(object.list = multiome_seurat_list)

# Integrate the datasets
multiome_seurat_integrated <- IntegrateData(anchorset = integration_anchors)

# Set the integrated assay as active
DefaultAssay(multiome_seurat_integrated) <- "integrated"

# Run standard analysis on integrated data
multiome_seurat_integrated <- ScaleData(multiome_seurat_integrated)
multiome_seurat_integrated <- RunPCA(multiome_seurat_integrated)
multiome_seurat_integrated <- RunUMAP(multiome_seurat_integrated, dims = 1:20)


int_author_midclust_plot = DimPlot(multiome_seurat_integrated, reduction = "umap", group.by = "mid_cluster", pt.size = 1) + 
  scale_color_manual(values = color_palette)
int_donor_plot = DimPlot(multiome_seurat_integrated, reduction = "umap", group.by = "orig.ident", pt.size = 1)

int_author_midclust_plot 
int_donor_plot 

DefaultAssay(multiome_seurat_integrated) <- "RNA"
gad2_p = FeaturePlot(multiome_seurat_integrated, features = "GAD2", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red")
gad1_p = FeaturePlot(multiome_seurat_integrated, features = "GAD1", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red")
vgat_p = FeaturePlot(multiome_seurat_integrated, features = "SLC32A1", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red")
vglut2_p = FeaturePlot(multiome_seurat_integrated, features = "SLC17A6", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red")
vglut1_p = FeaturePlot(multiome_seurat_integrated, features = "SLC17A7", reduction = "umap", pt.size = 1, slot = 'data') +
  scale_color_gradient(low = "white", high = "red")

#Look at the co-expression of specific genes
# Get expression data
umap_data <- as.data.frame(Embeddings(multiome_seurat_integrated, reduction = "umap"))
umap_data$GAD2 <- FetchData(multiome_seurat_integrated, vars = "GAD2", slot = "data")[, 1]
umap_data$SLC17A6 <- FetchData(multiome_seurat_integrated, vars = "SLC17A6", slot = "data")[, 1]
umap_data$SLC32A1 <- FetchData(multiome_seurat_integrated, vars = "SLC32A1", slot = "data")[, 1]

# Create a coexpression category
umap_data$coexpression <- ifelse(umap_data$GAD2 > 0 & umap_data$SLC17A6 > 0, "Both",
                                  ifelse(umap_data$GAD2 > 0, "GAD2 only",
                                         ifelse(umap_data$SLC17A6 > 0, "SLC17A6 only", "Neither")))

gad2_vglut2_p = ggplot(umap_data, aes(x = umap_1, y = umap_2, color = coexpression)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Both" = "purple", "GAD2 only" = "red", "SLC17A6 only" = "blue", "Neither" = "lightgrey")) +
  theme_bw() +
  labs(title = "GAD2 and SLC17A6 Co-expression")

gad2_p
gad1_p
vglut1_p
vglut2_p
gad2_vglut2_p

ggsave(plot = int_author_midclust_plot , path = plot_path, filename = 'LHb4_LHb7_integrated_midcluster_annot_umpa.pdf',
device = 'pdf', height = 5, width = 6, useDingbats = FALSE )

ggsave(plot = int_donor_plot , path = plot_path, filename = 'LHb4_LHb7_integrated_donor_annot_umpa.pdf',
device = 'pdf', height = 5, width = 6, useDingbats = FALSE )

ggsave(plot = gad1_p , path = plot_path, filename = 'LHb4_LHb7_integrated_gad1_annot_umap.pdf',
device = 'pdf', height = 5, width = 6, useDingbats = FALSE )


ggsave( plot = gad2_p, path = plot_path, filename = "LHb4_LHb7_integrated_GAD2_feature_umap.pdf",
  device = "pdf", height = 5, width = 6, useDingbats = FALSE
)

ggsave( plot = vgat_p, path = plot_path, filename = "LHb4_LHb7_integrated_SLC32A1_feature_umap.pdf",
  device = "pdf", height = 5, width = 6, useDingbats = FALSE
)

ggsave( plot = vglut1_p, path = plot_path, filename = "LHb4_LHb7_integrated_SLC17A7_feature_umap.pdf",
  device = "pdf", height = 5, width = 6, useDingbats = FALSE
)

ggsave( plot = vglut2_p, path = plot_path, filename = "LHb4_LHb7_integrated_SLC17A6_feature_umap.pdf",
  device = "pdf", height = 5, width = 6, useDingbats = FALSE
)

ggsave( plot = gad2_vglut2_p, path = plot_path, filename = "LHb4_LHb7_integrated_GAD2_SLC17A6_feature_umap.pdf",
  device = "pdf", height = 5, width = 6, useDingbats = FALSE
)

#Label the putative inhibitory LHb cells
DefaultAssay(multiome_seurat_integrated) <- "integrated"
# Find variable features
multiome_seurat_integrated <- FindVariableFeatures(multiome_seurat_integrated, selection.method = "vst", nfeatures = 2000)


# Find neighbors and clusters
multiome_seurat_integrated <- FindNeighbors(multiome_seurat_integrated, dims = 1:20)
multiome_seurat_integrated <- FindClusters(multiome_seurat_integrated, resolution = .3)

int_clust_p = DimPlot(multiome_seurat_integrated, reduction = "umap", group.by = "integrated_snn_res.0.3", pt.size = 1, label = TRUE)
int_clust_p

ggsave( plot =int_clust_p, path = plot_path, filename = "LHb4_LHb7_integrated_initial_new_cluster_umap.pdf",
  device = "pdf", height = 5, width = 6, useDingbats = FALSE
)

DefaultAssay(multiome_seurat_integrated) <- "RNA"
p_bubble = get_bubble_plot(multiome_seurat_integrated, 
  top_markers = c('CHAT', 'SLC5A7','SLC18A3', 'TAC1', 'TACR1','TAC3', 'GPR151', 'GAP43','SNAP25', 'POU4F1', 
  'SLC17A6', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1'),
 sample_name = "Multiome Habenula LHb4 and 7", group_col = "integrated_snn_res.0.3")

p_bubble[[1]]
p_bubble[[2]]

pdf(paste0(plot_path, '/Lhb4_7_integrated_multiome_bubble_Hab_marker_mean.pdf'), width = 10, height = 8)
p_bubble[[1]]
dev.off()

pdf(paste0(plot_path, '/Lhb4_7_integrated_multiome_bubble_Hab_marker_zscore_mean.pdf'), width = 10, height = 8)
p_bubble[[2]]
dev.off()

#Clearly cluster 5 snd 7

multiome_seurat_integrated$refined_mid_cluster = multiome_seurat_integrated$mid_cluster
multiome_seurat_integrated$refined_mid_cluster[multiome_seurat_integrated$integrated_snn_res.0.3 %in% c(5)] = 'Putative_Inhib_LHb_1'
multiome_seurat_integrated$refined_mid_cluster[multiome_seurat_integrated$integrated_snn_res.0.3 %in% c(7)] = 'Putative_Inhib_LHb_2'

table(multiome_seurat_integrated$refined_mid_cluster)
#LHb.4: 11704 cells
#LHb.7: 211 cells
#Putative_Inhib_LHb_1: 1126 cells
#Putative_Inhib_LHb_2: 871 cells


#Save the integrated seurat object

saveRDS(multiome_seurat_integrated, paste0(new_data_path, '/multiome_LHb4_LHb7_integrated_seurat.rds'))

inhib_meta = multiome_seurat_integrated[[]]


#Exclude them from the rest of the data set

multiome_sce




#Check out the top markers for the remaining cells, and then look at the expression of those markers 
#In the putative inhibitory cells.




#All the cell names are unique, so I can break apart by donor and then put it back together
table(duplicated(rownames(colData(multiome_sce))))

# Get unique donor IDs
donors <- unique(multiome_sce$orig.ident)

# Create a named list of SCE objects, one per donor
multiome_sce_list <- lapply(donors, function(donor) {
  multiome_sce[, multiome_sce$orig.ident == donor]
})
names(multiome_sce_list) <- donors

#Get the markers per donor
donor_markers_list = lapply(donors, function(donor) {
  MetaMarkers::compute_markers(assay(multiome_sce_list[[donor]], 'cpm'), multiome_sce_list[[donor]]$mid_cluster)
})
names(donor_markers_list) = donors

#Make metamarkers, really just interested in the thalamus markers
cross_donor_hab_markers = make_meta_markers(donor_markers_list, detailed_stats = TRUE)



clust_filt = multiome_sce[ , multiome_sce$mid_cluster == 'LHb.4']
thal_filt = multiome_sce[ , multiome_sce$mid_cluster == 'Excit.Thal']
in_thal_filt = multiome_sce[ , multiome_sce$mid_cluster == 'Inhib.Thal']

#This gets the cells that co-express GAD2 and SLC17A6
gad1_expression = assay(clust_filt, 'cpm')['GAD1', ]
gad2_expression = assay(clust_filt, 'cpm')['GAD2', ]
vglut2_expression = assay(clust_filt, 'cpm')['SLC17A6', ]
vglut1_expression = assay(clust_filt, 'cpm')['SLC17A7', ]

coGabaGlut_index = gad2_expression * vglut2_expression > 0
mean(coGabaGlut_index)
sum(coGabaGlut_index)
#Get the average expression per cell of the top thalamus markers
thalamus_markers = cross_donor_hab_markers %>% filter(cell_type %in% c('Excit.Thal', 'Inhib.Thal') & rank <= 25) %>% pull(gene)
hab_markers = cross_donor_hab_markers %>% filter(cell_type %in% c('LHb.4') & rank <= 50) %>% pull(gene)

avg_thal_expr = colMeans(assay(clust_filt, 'cpm')[thalamus_markers, ])
avg_hab_expr = colMeans(assay(clust_filt, 'cpm')[hab_markers, ])
gad2_exp = assay(clust_filt, 'cpm')['GAD2', ]
gad1_exp = assay(clust_filt, 'cpm')['GAD1', ]
coGabaGlut_label = rep('Not', length = length(coGabaGlut_index))
coGabaGlut_label[coGabaGlut_index] = 'CoGABA-Glut'
tissue_label = rep('Habenula', length = length(coGabaGlut_index))
co_exp_df_1 = data.frame(avg_thal_expr = avg_thal_expr, avg_hab_expr = avg_hab_expr, gad2_exp = gad2_exp, 
  coGabaGlut_label = coGabaGlut_label, tissue = tissue_label)

avg_thal_expr = colMeans(assay(thal_filt, 'cpm')[thalamus_markers, ])
avg_hab_expr = colMeans(assay(thal_filt, 'cpm')[hab_markers, ])
gad2_exp = assay(thal_filt, 'cpm')['GAD2', ]
gad1_exp = assay(thal_filt, 'cpm')['GAD1', ]
coGabaGlut_label = rep('Excite.Thal', length = ncol(thal_filt))
tissue_label = rep('Thalamus', length = ncol(thal_filt))
co_exp_df_2 = data.frame(avg_thal_expr = avg_thal_expr, avg_hab_expr = avg_hab_expr, gad2_exp = gad2_exp,
  coGabaGlut_label = coGabaGlut_label, tissue = tissue_label)

avg_thal_expr = colMeans(assay(in_thal_filt, 'cpm')[thalamus_markers, ])
avg_hab_expr = colMeans(assay(in_thal_filt, 'cpm')[hab_markers, ])
gad2_exp = assay(in_thal_filt, 'cpm')['GAD2', ]
gad1_exp = assay(in_thal_filt, 'cpm')['GAD1', ]
coGabaGlut_label = rep('Inhib.Thal', length = ncol(in_thal_filt))
tissue_label = rep('Thalamus', length = ncol(in_thal_filt))
co_exp_df_3 = data.frame(avg_thal_expr = avg_thal_expr, avg_hab_expr = avg_hab_expr, gad2_exp = gad2_exp,
  coGabaGlut_label = coGabaGlut_label, tissue = tissue_label)

co_exp_df = rbind(co_exp_df_1, co_exp_df_2,co_exp_df_3)
co_exp_df$coGabaGlut_label = factor(co_exp_df$coGabaGlut_label, levels = c('Not', 'CoGABA-Glut', 'Excite.Thal', 'Inhib.Thal'))
p_thal = ggplot(co_exp_df, aes(x = coGabaGlut_label, y = avg_thal_expr, fill = tissue)) +
  geom_violin(scale = 'width') +
  theme_bw() +
  xlab("Co-expression of GABA and Glut markers") +
  ylab("Average expression of top thalamus markers") +
  ggtitle("Thalamus marker expression in co-GABA-Glut cells vs others in LHb.4")

p_hab = ggplot(co_exp_df, aes(x = coGabaGlut_label, y = avg_hab_expr, fill = tissue)) +
  geom_violin(scale = 'width') +
  theme_bw() +
  xlab("Co-expression of GABA and Glut markers") +
  ylab("Average expression of top LHb.4 markers") +
  ggtitle("Habenula LHb.4 marker expression in co-GABA-Glut cells vs others in LHb.4")

p_gad = ggplot(co_exp_df, aes(x = coGabaGlut_label, y = gad2_exp, fill = tissue)) +
  geom_violin(scale = 'width') +
  theme_bw() +
  xlab("Co-expression of GABA and Glut markers") +
  ylab("GAD2 expression") +
  ggtitle("GAD2 expression in co-GABA-Glut cells vs others in LHb.4")

p_thal
p_hab
p_gad


############
#Okay so not obviously thalamus cells based on thalamus marker expression
###########

ggsave(plot = p_thal, filename = 'Thalamus_marker_expression_by_GABA_Glut_coexpression_status.pdf', path = plot_path, 
device = 'pdf', width = 6, height = 5, useDingbats = FALSE)

ggsave(plot = p_hab, filename = 'habenula_lhb4_marker_expression_by_GABA_Glut_coexpression_status.pdf', path = plot_path, 
device = 'pdf', width = 6, height = 5, useDingbats = FALSE)

ggsave(plot = p_gad, filename = 'GAD2_expression_by_GABA_Glut_coexpression_status.pdf', path = plot_path, 
device = 'pdf', width = 6, height = 5, useDingbats = FALSE)

