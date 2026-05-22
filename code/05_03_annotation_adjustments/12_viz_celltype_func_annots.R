#exploring and plotting the gene set annotation results

library(MetaNeighbor)
library(MetaMarkers)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(qs2)
library(here)
library(sessioninfo)

message('Loading data...')

multiome_path = here('processed-data', '05_03_annotation_adjustments','06_refined_annotations')
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

assay(multiome_sce, 'cpm') = MetaMarkers::convert_to_cpm(assay(multiome_sce, 'counts'))

#Filter for just the habenula cell-types
hab_celltypes = c('Inhib_LHb_4.1','Inhib_LHb_4.2','LHb.4','LHb.2.7','LHb.1.3.4','MHb.1','MHb.1.2','MHb.2','MHb.3')
multiome_sce = multiome_sce[ ,multiome_sce$refined_mid_cluster %in% hab_celltypes]


#Drop ATAC assay
altExps(multiome_sce) <- NULL
gc()

kegg_path = here('processed-data','05_03_annotation_adjustments', '11_kegg_term_annot')
go_path = here('processed-data','05_03_annotation_adjustments', '10_celltype_func_annot')

#GO and KEGG results
go_aurocs = read.table(here(go_path, 'multiomeHab_functional_aurocs.txt'))
kegg_aurocs = read.table(here(kegg_path, 'multiomeHab_KEGG_functional_aurocs.txt'))

#GO and KEGG terms
kegg_sets = readRDS(file.path(kegg_path, "kegg_human.rds"))
go_sets = readRDS(file.path(go_path, "go_human.rds"))

plot_path = here('plots','05_03_annotation_adjustments', '12_viz_celltype_func_annots')
if (!dir.exists(plot_path)) dir.create(plot_path)


#Set up kegg set stats
gs_size = sapply(kegg_sets, length)
kegg_aurocs_df = data.frame(kegg_term = rownames(kegg_aurocs), kegg_aurocs)
kegg_aurocs_df$average = rowMeans(kegg_aurocs)
kegg_aurocs_df$n_genes = gs_size[rownames(kegg_aurocs)]
kegg_aurocs_df$ontology = 'KEGG'

kegg_aurocs_df


head(kegg_aurocs_df[order(kegg_aurocs_df$average, decreasing = TRUE),], 10)

small_kegg_sets = kegg_aurocs_df[kegg_aurocs_df$n_genes <= 50,]
head(small_kegg_sets[order(small_kegg_sets$average, decreasing = TRUE),], 20) 

set_of_interest = 'Aldosterone-regulated sodium reabsorption'
plotDotPlot(dat = multiome_sce,
experiment_labels = multiome_sce$orig.ident,
celltype_labels = multiome_sce$refined_mid_cluster,
gene_set = kegg_sets[[set_of_interest]]) + ggtitle(set_of_interest)


#Highlight the high scores of the neuropsychiatric terms
terms_to_label = c('Amphetamine addiction','Long-term depression','Nicotine addiction','Morphine addiction','Cocaine addiction')
kegg_aurocs_df$label_terms = ifelse(kegg_aurocs_df$kegg_term %in% terms_to_label, kegg_aurocs_df$kegg_term, NA)
kegg_aurocs_df$point_color = ifelse(kegg_aurocs_df$kegg_term %in% terms_to_label, 'red', 'black')
pos <- position_jitter(width = 0.1, seed = 2)

ggplot(kegg_aurocs_df, aes(x = ontology, y = average, label = label_terms)) + 
  geom_boxplot(fill = 'grey') + 
  ggrepel::geom_label_repel(min.segment.length = 0, max.overlaps = Inf,
    box.padding = 1, position = pos,
    xlim  = c(NA ,1) ) +
  geom_point(position = pos, aes(color = point_color), show.legend = FALSE) + 
  scale_color_manual(values = c('red' = 'red', 'black' = 'black')) +
  theme_bw() + ylab('Predictive strength (avg. AUROC)') +
  ggtitle('KEGG gene sets predicting Habenula cell-types')






#Set up go set stats
gs_size = sapply(go_sets, length)
go_aurocs_df = data.frame(go_term = rownames(go_aurocs), go_aurocs)
go_aurocs_df$average = rowMeans(go_aurocs)
go_aurocs_df$n_genes = gs_size[rownames(go_aurocs)]
go_aurocs_df$ontology = 'GO'
head(go_aurocs_df[order(go_aurocs_df$MHb.1, decreasing = TRUE),], 20)

small_go_sets = go_aurocs_df[go_aurocs_df$n_genes <= 30,]
head(small_go_sets[order(small_go_sets$Inhib_LHb_4.2, decreasing = TRUE),], 20) 

set_of_interest = "GO:1903859|regulation of dendrite extension|BP"
plotDotPlot(dat = multiome_sce,
experiment_labels = multiome_sce$orig.ident,
celltype_labels = multiome_sce$refined_mid_cluster,
gene_set = go_sets[[set_of_interest]]) + ggtitle(set_of_interest)



#Highlight the high scores of the neuropsychiatric terms
terms_to_label = go_aurocs_df %>% slice_max(order_by = average, n = 10) %>% pull(go_term)
go_aurocs_df$label_terms = ifelse(go_aurocs_df$go_term %in% terms_to_label, go_aurocs_df$go_term, NA)
go_aurocs_df$point_color = ifelse(go_aurocs_df$go_term %in% terms_to_label, 'red', 'black')
pos <- position_jitter(width = 0.1, seed = 2)

ggplot(go_aurocs_df, aes(x = ontology, y = average, label = label_terms)) + 
  geom_boxplot(fill = 'grey') + 
  ggrepel::geom_label_repel(min.segment.length = 0, max.overlaps = Inf,
    box.padding = 1, position = pos,
    xlim  = c(NA ,1) ) +
  geom_point(position = pos, aes(color = point_color), show.legend = FALSE) + 
  scale_color_manual(values = c('red' = 'red', 'black' = 'black')) +
  theme_bw() + ylab('Predictive strength (avg. AUROC)') +
  ggtitle('go gene sets predicting Habenula cell-types')




