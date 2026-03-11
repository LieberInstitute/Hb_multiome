#The LHb.4 lateral population is the one that maps to the mouse and zebrafish lateral clusters that express both GABA and Glut markers
#Diving into that feature a bit more


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
new_data_path = here('processed-data', '99_donor_cluster_replicability','01_per_donor_clustering')
#Path to plot directory
plot_path = here('plots', '99_donor_cluster_replicability','01_per_donor_clustering')

if (!dir.exists(new_data_path)) dir.create(new_data_path)
if (!dir.exists(plot_path)) dir.create(plot_path)


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


############################
#While the LHb.4 cluster does have some co-expression of GABA and Glut markers, they are super small percentages.
#I'm curious if these are just small subsets of thalamus cells that are present within the LHb.4 cluster
#Two pieces of evidence, the LHb.4 cluster is spatially the one closest to the thalamus
#And the GABA glut co-expression is the strongest in the thalamus clusters.

#So even a small amount of thalamus in the LHb.4 could account for the low levels of co GABA glut

##############################

#For the subset of GABA-Glut cells in a cluster, compare the expression of thalamus markers for those cells

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
ggplot(co_exp_df, aes(x = coGabaGlut_label, y = avg_thal_expr, fill = tissue)) +
  geom_violin(scale = 'width') +
  theme_bw() +
  xlab("Co-expression of GABA and Glut markers") +
  ylab("Average expression of top thalamus markers") +
  ggtitle("Thalamus marker expression in co-GABA-Glut cells vs others in LHb.4")

ggplot(co_exp_df, aes(x = coGabaGlut_label, y = avg_hab_expr, fill = tissue)) +
  geom_violin(scale = 'width') +
  theme_bw() +
  xlab("Co-expression of GABA and Glut markers") +
  ylab("Average expression of top LHb.4 markers") +
  ggtitle("Thalamus marker expression in co-GABA-Glut cells vs others in LHb.4")

ggplot(co_exp_df, aes(x = coGabaGlut_label, y = gad2_exp, fill = tissue)) +
  geom_violin(scale = 'width') +
  theme_bw() +
  xlab("Co-expression of GABA and Glut markers") +
  ylab("GAD2 expression") +
  ggtitle("GAD2 expression in co-GABA-Glut cells vs others in LHb.4")



############
#Okay so not obviously thalamus cells based on thalamus marker expression
###########